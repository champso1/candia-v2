#include "Candia-v2/Grid.hpp"
#include "Candia-v2/Common.hpp"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_interp.h>
#include <iterator>
#include <limits>
#include <set>
#include <ranges>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace Candia2
{
	Grid::Grid(grid_type const& xtab, uint nx, uint grid_fill_type)
		: _points(nx), _ntab{}, _xtab{xtab}, _interval_sizes{DEFAULT_GAULEG_POINTS},
		  _Xi(1, gauleg_type(DEFAULT_GAULEG_POINTS, 0.0)), _Wi(1, gauleg_type(DEFAULT_GAULEG_POINTS, 0.0))
	{
		if (!std::ranges::is_sorted(xtab)) {
			log(LOG_WARNING, "Grid", "xtab array was not sorted. will sort it and continue");
			std::ranges::sort(_xtab);
		}

		_workspaces.reserve(DISTS);
		_interps.reserve(DISTS);
		_interp_accels.reserve(DISTS);
		for (uint i=0; i<DISTS; ++i)
			_workspaces.emplace_back(gsl::make_default_workspace());

	    assert(gsl::error_handler_set, "GSL error handler failed to set correctly.");
		
		switch (grid_fill_type)
		{
			case LOG: initGridLog(_xtab, nx); break;
			case LOG_LIN: initGridLogLin(_xtab, nx); break;
			case LIN: initGridLin(_xtab, nx); break;
			case LOG_LIN_QUAD: initGridLogLinQuad(_xtab, nx); break;
			default:
				log(LOG_ERROR, "Grid", "Invalid grid fill type. Found {}, expected 1(LOG) or 2(LOG_LIN).", grid_fill_type);
		}
		
		initGauLeg(0.0, 1.0, _Xi[0], _Wi[0]);
	}

	Grid::Grid(Grid& other)
		: OptionsBase<GridOptions>{other},
		_points{other._points}, _ntab{other._ntab}, _xtab{other._xtab},
		_split_interval{other._split_interval},
		_Xi{other._Xi}, _Wi{other._Wi},
		_interval_sizes{other._interval_sizes}
	{
		_workspaces.reserve(DISTS);
		_interps.reserve(DISTS);
		_interp_accels.reserve(DISTS);
		for (uint i=0; i<DISTS; ++i)
			_workspaces.emplace_back(gsl::make_workspace(other._workspaces[i]->limit));

		if (!other._interps.empty()) {
			for (uint i=0; i<DISTS; ++i) {
				_interps.emplace_back(gsl::make_interp(other._interps[i]->size));
				_interp_accels.emplace_back(gsl::make_interp_accel());
			}
		}
	}
	
	void Grid::operator=(Grid& other)
	{
		_points = other._points;
		_ntab = other._ntab;
		_xtab = other._xtab;
		_split_interval = other._split_interval;
		_Xi = other._Xi;
		_Wi = other._Wi;
		_interval_sizes = other._interval_sizes;
		options = other.options;

		_workspaces.clear();
		_interps.clear();
		_interp_accels.clear();
		
		_workspaces.reserve(DISTS);
		_interps.reserve(DISTS);
		_interp_accels.reserve(DISTS);
		for (uint i=0; i<DISTS; ++i)
			_workspaces.emplace_back(gsl::make_workspace(other._workspaces[i]->limit));

		if (!other._interps.empty()) {
			for (uint i=0; i<DISTS; ++i) {
				_interps.emplace_back(gsl::make_interp(other._interps[i]->size));
				_interp_accels.emplace_back(gsl::make_interp_accel());
			}
		}
	}

	void Grid::initGridLog(grid_type const& xtab, uint nx)
	{
		const uint xtab_len = xtab.size();
		std::vector<double> Ntab(xtab_len);
	    std::vector<int> ntab(xtab_len);

		double temp = -std::log10(xtab[0]);

		for (uint i=1; i<xtab_len; i++)
			Ntab[i] = (double)(nx-1)*std::log10(xtab[i]/xtab[i-1])/temp;

		ntab[0] = nx-1;

		for (uint i=1; i<xtab_len; i++) {
			ntab[i] =  (int)Ntab[i];
			Ntab[i] -= (double)ntab[i];
			ntab[0] -= ntab[i];
		}

		for (uint i=1; i<xtab_len; i++) {
			if (ntab[i] == 0) {
				ntab[i] =  1;
				Ntab[i] -= 1.0;
				ntab[0] -= 1;
			}
		}

		uint n;
		for ( ; ntab[0]<0; ntab[0]++) {
			n=0;
			for (uint i=1; i<xtab_len; i++) {
				if (ntab[i] != 1) {
					if ((n == 0) || (Ntab[i] <= temp)) {
						n = i;
						temp = Ntab[i];
					}
				}
			}

			ntab[n]--;
			Ntab[n] += 1.0;
		}

		for ( ; ntab[0]>0; ntab[0]--) {
			n=0;
			for (uint i=1; i<xtab_len; i++) {
				if ((n == 0) || (Ntab[i] > temp)) {
					n = i;
					temp = Ntab[i];
				}
			}

			ntab[n]++;
			Ntab[n] -= 1.0;
		}

		for (uint i=1; i<xtab_len; i++)
			ntab[i] += ntab[i-1];

		double lstep;
		for (uint i=0; i<xtab_len-1; i++) {
			lstep = std::log10(xtab[i+1]/xtab[i])/(double)(ntab[i+1]-ntab[i]);

			for (int j=ntab[i]; j<ntab[i+1]; j++)
				_points.at(j) = xtab[i]*std::pow(10.0, lstep*(double)(j-ntab[i]));
		}

		_points.at(nx-1) = 1.0;
		_ntab = ntab;

		return;
	}



	void Grid::initGridLogLin(grid_type const& xtab, uint nx)
	{
		// this is to make the intervals less "clean"
		// sometimes, when they are "clean", the linear mapping places points basically right on
		// the xtab points, like 0.3, but off by a delta which is small enough to mess up interpolation
		// if this number is a bit uneven, the hope is that points won't be placed so "cleanly" near
		// xtabbed points, avoiding these errors
		nx += 1;

		double log_min = std::log10(1e-5);
		double log_max = std::log10(0.1);
		uint num_log_intervals = std::round(log_max-log_min);
		double dlog = (log_max-log_min)/static_cast<double>(num_log_intervals);
		uint log_interval_size = nx/num_log_intervals;

		double lin_min = 0.1;
		double lin_max = 1.0;

		_points.clear();
		for (uint i=0; i<num_log_intervals; ++i) {
			double l0 = log_min + i*dlog;
			double l1 = l0 + dlog;
			for (uint k=0; k<log_interval_size; ++k) {
				double l = l0 + (l1-l0)*k/static_cast<double>(log_interval_size);
				_points.emplace_back(std::pow(10, l));
			}
		}
		
		for (uint k=0; k<nx; ++k) {
		    double x = lin_min + (lin_max-lin_min)*k/static_cast<double>(nx);
			_points.emplace_back(x);
		}
		
		std::set<double> points_set(_points.begin(), _points.end());
		points_set.insert(xtab.begin(), xtab.end());
		_points = std::vector<double>(points_set.begin(), points_set.end());
		
		_ntab.clear();
		for (double x : xtab) {
			if (auto it = std::find(_points.begin(), _points.end(), x); it != _points.end()) {
				_ntab.emplace_back(std::distance(_points.begin(), it));
				continue;
			}
		}
	}

	void Grid::initGridLin(grid_type const& xtab, uint nx)
	{
		// this is to make the intervals less "clean"
		// sometimes, when they are "clean", the linear mapping places points basically right on
		// the xtab points, like 0.3, but off by a delta which is small enough to mess up interpolation
		// if this number is a bit uneven, the hope is that points won't be placed so "cleanly" near
		// xtabbed points, avoiding these errors
		nx += 1;
		double lin_min = 1e-5;
		double lin_max = 1.0;

		_points.clear();
		
		for (uint k=0; k<nx; ++k) {
		    double x = lin_min + (lin_max-lin_min)*k/static_cast<double>(nx);
			_points.emplace_back(x);
		}
		
		std::set<double> points_set(_points.begin(), _points.end());
		points_set.insert(xtab.begin(), xtab.end());
		_points = std::vector<double>(points_set.begin(), points_set.end());
		
		_ntab.clear();
		for (double x : xtab) {
			if (auto it = std::find(_points.begin(), _points.end(), x); it != _points.end()) {
				_ntab.emplace_back(std::distance(_points.begin(), it));
				continue;
			}
		}
	}

	void Grid::initGridLogLinQuad(grid_type const& xtab, uint nx)
	{
		// this is to make the intervals less "clean"
		// sometimes, when they are "clean", the linear mapping places points basically right on
		// the xtab points, like 0.3, but off by a delta which is small enough to mess up interpolation
		// if this number is a bit uneven, the hope is that points won't be placed so "cleanly" near
		// xtabbed points, avoiding these errors
		nx += 1;

		double log_min = std::log10(1e-5);
		double log_max = std::log10(0.1);
		uint num_log = 51;
		uint num_log_intervals = std::round(log_max-log_min);
		double dlog = (log_max-log_min)/static_cast<double>(num_log_intervals);
		uint log_interval_size = num_log/num_log_intervals;

		_points.clear();
		for (uint i=0; i<num_log_intervals; ++i) {
			double l0 = log_min + i*dlog;
			double l1 = l0 + dlog;
			for (uint k=0; k<log_interval_size; ++k) {
				double l = l0 + (l1-l0)*k/static_cast<double>(log_interval_size);
				_points.emplace_back(std::pow(10, l));
			}
		}

		double lin_min = 0.1;
		double lin_max = 0.9;
		uint num_lin = 26;
		
		for (uint k=0; k<num_lin; ++k) {
		    double x = lin_min + (lin_max-lin_min)*k/static_cast<double>(num_lin);
			_points.emplace_back(x);
		}

		double quad_min = 0.9;
		double quad_max = 1.0;
		uint num_quad = 26;

		for (uint k=0; k<num_quad; ++k) {
			double f = 1.0 - k/static_cast<double>(num_quad);
		    double x = quad_max - (quad_max-quad_min)*f*f;
			_points.emplace_back(x);
		}
		
		std::set<double> points_set(_points.begin(), _points.end());
		points_set.insert(xtab.begin(), xtab.end());
		_points = std::vector<double>(points_set.begin(), points_set.end());
		
		_ntab.clear();
		for (double x : xtab) {
			if (auto it = std::find(_points.begin(), _points.end(), x); it != _points.end()) {
				_ntab.emplace_back(std::distance(_points.begin(), it));
				continue;
			}
		}

		for (uint i=0; i<DISTS; ++i) {
			_interps.emplace_back(gsl::make_interp(_points.size()));
			_interp_accels.emplace_back(gsl::make_interp_accel());
		}
	}

	void Grid::initGauLeg(value_type x1, value_type x2, gauleg_type& Xi, gauleg_type& Wi)
	{
		const double eps = 3.0e-11; // relative precision

		// abscissae are symmetric:
		uint n = Xi.size(); // simpler to type
		double N = static_cast<double>(n);
		uint m = (n+1)/2;
		// double x2 = 1.0;
		// double x1 = 0.0;
		double xm = 0.5*(x2+x1);
		double xl = 0.5*(x2-x1);

		for (uint i=1; i<=m; i++) {
			double I = static_cast<double>(i);
			double z = std::cos(PI*(I-0.25)/(N+0.5));

			// default initialize some of these.
			// easy to spot if error occurs
			double z1 = std::numeric_limits<double>::max();
			double pp = std::numeric_limits<double>::max();

			double p1, p2, p3;
			double J;
			do {
				p1 = 1.0;
				p2 = 0.0;

				for (uint j=1; j<=n; j++) {
					J = static_cast<double>(j);
					p3 = p2;
					p2 = p1;
					p1 = ((2.0*J - 1.0)*z*p2 - (J-1.0)*p3)/J;
				}

				pp = N*(z*p1 - p2)/(z*z - 1.0);
				z1 = z;
				z = z1 - p1/pp;
			} while (std::abs(z-z1) > eps);

			if (z1 == std::numeric_limits<double>::max() ||
				pp == std::numeric_limits<double>::max())
				log(LOG_ERROR, "Grid::initGauLeg()", "Failed to determine gauss-legendre abscissae/weights");

			Xi[i-1] = xm - xl*z;
			Xi[n-i] = xm + xl*z;
			Wi[i-1] = 2.0*xl/((1.0 - z*z)*pp*pp);
			Wi[n-i] = Wi[i-1];
		}
	}

	void Grid::splitConvolution(
		std::vector<value_type> const& intervals,
		std::vector<value_type> const& sizes)
	{
		// this means we just accept the default splitting
		if (intervals.empty() && sizes.empty()) {
			log(LOG_INFO, "Grid", "Using default set of intervals");
			_split_interval = true;
		    return;
		} else if ((intervals.empty() && !sizes.empty()) || (!intervals.empty() && sizes.empty()))
			log(LOG_ERROR, "Grid::splitConvolution()", "Must provide either both intervals and sizes or neither (for default).");
		
		if (intervals.size()-1 != sizes.size())
			log(LOG_ERROR, "Grid::splitConvolution()", "Invalid number of element in the intervals and sizes vectors.");
		
		auto num_intervals = intervals.size()-1;
		if (num_intervals > 5)
			log(LOG_WARNING, "Grid::splitConvolution()", ">5 convolution intervals is a little excessive.");
		
		if (auto it = std::find_if(intervals.begin(), intervals.end(), [](double x) { return (x < 1e-7) || (x > 1.0); });
			it != std::ranges::end(intervals))
		{
			log(LOG_ERROR, "Grid::splitConvolution()", "A provided point ({}) is outside the range [1e-7, 1.0]", *it);
		}

		_Xi.clear();
		_Wi.clear();
		_interval_sizes.clear();
		{
			gauleg_type X(DEFAULT_GAULEG_POINTS, 0.0), W(DEFAULT_GAULEG_POINTS, 0.0);
			initGauLeg(0, 1.0, X, W);
			_Xi.emplace_back(X);
			_Wi.emplace_back(W);
		}
		
		auto points_per_interval = size()/num_intervals;
		
		auto it = intervals.begin() + 1;
		auto it_size = sizes.begin();
		while (it != intervals.end() && it_size != sizes.end()) {
			auto prev_x = *(it-1);
			auto new_x = *it;
			auto size = *it_size;
			gauleg_type Xi(size, 0.0), Wi(size, 0.0);
			initGauLeg(prev_x, new_x, Xi, Wi);
			_Xi.emplace_back(Xi);
			_Wi.emplace_back(Wi);
			_interval_sizes.emplace_back(size);
			
			++it;
			++it_size;
		}

		log(LOG_DEBUG, "Grid::splitConvolution", "Using the following sizes: {}", vec_to_str(_interval_sizes));
		auto tie =
			std::views::iota(0)
			| std::views::take(num_intervals)
			| std::views::transform([&](int i){ return std::make_tuple(i, _Xi[i+1], _Wi[i+1], _interval_sizes[i]); });
		for (auto const& [i, X, W, s] : tie) {
			auto min = intervals[i];
			auto max = intervals[i+1];
			
			log(LOG_DEBUG, "Grid::splitConvolution()", "Interval of size {} in range [{}, {}]", s, min, max);
		}
		
		assert(_Xi.size() == _Wi.size() && _Xi.size() == _interval_sizes.size()+1, "Failed to split convolution intervals.");
		_split_interval = true;
	}

    int Grid::interpFindIdx(value_type x)
	{
		const static int n = static_cast<int>(size());
		const static int max_k = static_cast<int>(n-2*INTERP_POINTS);

		auto it = std::upper_bound(_points.begin(), _points.end(), x);
		int k = static_cast<int>(it - _points.begin()) - INTERP_POINTS;
		
		if (k < 0)
			k = 0;
		else if (k > max_k)
			k = max_k;

		return k;
	}

	Grid::value_type Grid::interpolate(arraygrid_type& yy, value_type x)
	{
		const static int n = 2*INTERP_POINTS;
		int ns=0;
		double y, den, dif, dift, ho, hp, w;

		int k = interpFindIdx(x);

		
		double const* xa = &(_points.data()[k]);
		double const* ya = &(yy.base().data()[k]);
		
		std::array<double, 2*INTERP_POINTS> c{};
		std::array<double, 2*INTERP_POINTS> d{};

		dif = std::abs(x - xa[0]);

		for (int i=0; i<n; i++) {
			if ((dift = std::abs(x - xa[i])) < dif) {
				ns = i;
				dif = dift;
			}
			c[i] = ya[i];
			d[i] = ya[i];
		}

		y = ya[ns--];

		for (int m=1; m<n; m++) {
		    for (int i=0; i<n-m; i++) {
				ho = xa[i] - x;
				hp = xa[i+m] - x;
				w = c[i+1] - d[i];

				den = ho-hp;
				if (std::abs(ho-hp) < 1e-15) {
					log(LOG_ERROR_NOQUIT, "Grid::interpolate()", "found a denominator equal to 0.0 ({}, {}, {})", x, xa[i], xa[i+m]);
					std::string msg = std::format("found a denominator equal to 0.0 ({}, {}, {})", x, xa[i], xa[i+m]);
					throw std::runtime_error("interpolate failed.");
				}

				den = w/den;
				d[i] = hp*den;
				c[i] = ho*den;
			}

			y += (2*ns < (n-1-m) ? c[ns+1] : d[ns--]);
		}

		return y;
	}

	static double gsl_convolution_function(double y, void* params_)
	{
		Grid::GSLIntegrationParams* params = reinterpret_cast<Grid::GSLIntegrationParams*>(params_);
		auto& g = params->g;
		auto x = params->x;
		auto k = params->k;
		auto eplus1 = params->eplus1;
		auto& A = params->A;
		auto& E = params->E;

		double Ak = A[k];
		double a = x/y;
		
		double res = 0.0;
		{ // regular piece
			double interpy = g.interpolate(A, y);
			res += (1.0/y)*a*interpy;
		}
		{ // delta function piece
			double interpa = g.interpolate(A, a);
			double eplusy = E.calcPlus(y);
			res += (eplusy*interpa - eplus1*Ak)/(1-y);
		}
		
		return res;
	}

	Grid::value_type Grid::mappingFunctionBase(
		uint k, value_type x, auto&& yandjaccessor,
		expression_type& E, arraygrid_type& A,
		value_type eplus1,
		gauleg_type const& X, gauleg_type const& W)
	{
		// this gsl stuff must be reset everytime we perform a new convolution
		auto idx = thread_index+1;
		auto* interp = _interps[idx].get();
		auto* accel = _interp_accels[idx].get();
		auto const* xa = _points.data();
		auto const* ya = A.base().data();
		gsl_interp_accel_reset(accel);
		gsl_interp_init(interp, xa, ya, size());
		
		double ak = A[k];
		double out = 0.0;
		uint s = X.size();
		for (uint i=0; i<s; i++) {
			double z = X[i];
			double w = W[i];
			auto [y, J] = yandjaccessor(x, z);
			double a = x/y;
			double eregy = E.calcRegular(y);
			double interpa = gsl_interp_eval(interp, xa, ya, a, accel);
			// double interpa = interpolate(A, a);
			double eplusy = E.calcPlus(y);
			
			out += w*J * eregy*interpa;
			out += w*J * (1.0/(1.0-y))*(eplusy*interpa - eplus1*ak);
		}
		return out;
	}

	Grid::value_type Grid::convolution(arraygrid_type& A, expression_type &E, uint k)
	{
		static int print_count = 100;

		double x = _points[k];
		double d1 = 1.0-x;
		double d2 = x*d1;
		double logx =  std::log(x);
		double eplus1 = E.calcPlus(1.0);
		double ed1 = E.calcDelta(1.0);
		double res = (eplus1*std::log1p(-x) + ed1) * A[k];

		auto gauleg_conv = [&](gauleg_type const& X, gauleg_type const& W) {
			for (uint i=0; i<X.size(); i++) {
				double z = X[i];
				double w = W[i];

				double y = std::pow(x, z);
				double a = std::pow(x, 1.0-z);

				double interpa = interpolate(A, a);

				double eregy = E.calcRegular(y);
				double eplusy = E.calcPlus(y);

				res -= w*logx * y*eregy*interpa;
				res -= w*logx * y*(eplusy*interpa - eplus1*A[k])/(1.0-y);
			}
		};

		// we use this nice convenience view for skipping the first gauleg points/weights
		// which are the default 0-1 that is used directly
		static auto gauleg_enum = [&](uint d) {
			return 
				std::views::iota(1)
				| std::views::take(_Xi.size()-1-d)
				| std::views::transform([&](int i){ return std::make_tuple(i, _Xi.begin()+i, _Wi.begin()+i); });
		};
		
		if (!_split_interval) {
			gauleg_conv(_Xi[0], _Wi[0]);
		} else {
			if (!options.use_gsl_conv_routine && !options.use_alt_mapping) {
				for (auto const& [i, X, W] : gauleg_enum(0))
					gauleg_conv(*X, *W);

			} else if (!options.use_gsl_conv_routine && options.use_alt_mapping) {
				// [x, 0.1]
			    static auto mapping1 = [&](double x, double z){
					auto a = 0.1/x;
					return std::make_pair(x*std::pow(a, z), x*std::pow(a, z)*std::log(a)); };

				// [0.1,0.9]
				static auto mapping2 = [&](double x, double z){
					return std::make_pair(0.1+0.8*z, 0.8); };
				static auto mapping2x = [&](double x, double z){
					return std::make_pair(x+(0.9-x)*z, 0.9-x); };

				// [0.9, 1.0]
				static auto mapping3 = [&](double x, double z){
					return std::make_pair(1.0-0.1*std::pow(1.0-z, 3), 0.3*(1.0-z)*(1.0-z)); };
				static auto mapping3x = [&](double x, double z){
					return std::make_pair(1.0-(1.0-x)*std::pow(1.0-z, 3), 3.0*(1-x)*(1.0-z)*(1.0-z)); };

				if (x < 0.1) {
					res += mappingFunctionBase(k, x, mapping1, E, A, eplus1, _Xi[0], _Wi[0]);
					res += mappingFunctionBase(k, x, mapping2, E, A, eplus1, _Xi[0], _Wi[0]);
					res += mappingFunctionBase(k, x, mapping3, E, A, eplus1, _Xi[0], _Wi[0]);
				} else if (x >= 0.1 && x < 0.9) {
					res += mappingFunctionBase(k, x, mapping2x, E, A, eplus1, _Xi[0], _Wi[0]);
					res += mappingFunctionBase(k, x, mapping3, E, A, eplus1, _Xi[0], _Wi[0]);
				} else {
					if ((1.0-x) < 1e-10)
						res += 0;
					else
						res += mappingFunctionBase(k, x, mapping3x, E, A, eplus1, _Xi[0], _Wi[0]);
				}
			} else {
				GSLIntegrationParams p{
					.g = *this,
					.x = x,
					.k = k,
					.logx = logx,
					.eplus1 = eplus1,
					.A = A,
					.E = E};
				
				gsl_function f{
					.function = gsl_convolution_function,
					.params = reinterpret_cast<void*>(&p) };
				double a = x, b = 1.0-1e-15;
				double epsabs = 0.0, epsrel = 1e-5;
				int key = GSL_INTEG_GAUSS21;
				std::size_t limit = gsl::DEFAULT_WORKSPACE_SIZE;
				double out{}, abserr{};
				gsl_integration_workspace* w = _workspaces[thread_index+1].get();

				int rc = gsl_integration_qags(
					&f, a, b, epsabs, epsrel,
					limit, w,
					&out, &abserr);
				if (rc != GSL_SUCCESS) {
				    _gsl_conv_errors.emplace_back(x, out, abserr);
				}
				
				res += out;
			}
		}
		
		return res;
	}
}
