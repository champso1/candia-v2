#include "Candia-v2/Grid.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/ArrayGrid.hpp"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <iterator>
#include <limits>
#include <set>
#include <ranges>

namespace Candia2
{
	Grid::Grid(std::vector<double> const& xtab, uint nx, uint grid_fill_type)
		: _points(nx), _ntab{}, _xtab{xtab}, _interval_sizes{DEFAULT_GAULEG_POINTS},
		  _Xi{gauleg_type(DEFAULT_GAULEG_POINTS, 0.0)}, _Wi{1, gauleg_type(DEFAULT_GAULEG_POINTS, 0.0)},
		  _workspace{gsl::make_default_workspace()}
	{
		if (!std::ranges::is_sorted(xtab)) {
			log(LOG_WARNING, "Grid", "xtab array was not sorted. will sort it and continue");
			std::ranges::sort(_xtab);
		}

		_workspaces.reserve(DISTS);
		for (uint i=0; i<DISTS; ++i)
			_workspaces.emplace_back(gsl::make_default_workspace());

	    assert(gsl::error_handler_set, "GSL error handler failed to set correctly.");
		
		switch (grid_fill_type)
		{
			case LOG: initGridLog(_xtab, nx); break;
			case LOG_LIN: initGridLogLin(_xtab, nx); break;
			default:
				log(LOG_ERROR, "Grid", "Invalid grid fill type. Found {}, expected 1(LOG) or 2(LOG_LIN).", grid_fill_type);
		}
		
		initGauLeg(0.0, 1.0, _Xi[0], _Wi[0]);
	}

	Grid::Grid(Grid& other)
		: _points{other._points}, _ntab{other._ntab}, _xtab{other._xtab},
		  _split_interval{other._split_interval},
		  _Xi{other._Xi}, _Wi{other._Wi},
		  _interval_sizes{other._interval_sizes},
		  _use_gsl_routine_for_high_x(other._use_gsl_routine_for_high_x),
		  _workspace{std::move(gsl::make_workspace(other._workspace->limit))}
	{
		_workspaces.reserve(DISTS);
		for (uint i=0; i<DISTS; ++i)
			_workspaces.emplace_back(gsl::make_workspace(other._workspaces[i]->limit));
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
		_use_gsl_routine_for_high_x = other._use_gsl_routine_for_high_x;
	    _workspace = std::move(gsl::make_workspace(other._workspace->limit));

		_workspaces.clear();
		_workspaces.reserve(DISTS);
		for (uint i=0; i<DISTS; ++i)
			_workspaces.emplace_back(gsl::make_workspace(other._workspaces[i]->limit));
	}

	void Grid::initGridLog(std::vector<double> const& xtab, const uint nx)
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



	void Grid::initGridLogLin(std::vector<double> const& xtab, uint nx)
	{
		log(LOG_WARNING, "Grid", "Method 2 is unfinished. Prefer method 3 for now.");
		log(LOG_WARNING, "Grid", "Will ignore supplied x-tab. Supplied number of grid points will actually be 1/3 of total grid points.");

		double log_min = std::log10(1e-5);
		double log_max = std::log10(0.1);
		uint num_log_intervals = std::round(log_max-log_min);
		double dlog = (log_max-log_min)/static_cast<double>(num_log_intervals);
		uint log_interval_size = nx/num_log_intervals;
		
		double lin1_min = 0.1;
		double lin1_max = 0.8;
		double lin2_min = 0.8;
		double lin2_max = 1.0;

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
		    double x = lin1_min + (lin1_max-lin1_min)*k/static_cast<double>(nx);
			_points.emplace_back(x);
		}
		for (uint k=0; k<nx; ++k) {
		    double x = lin2_min + (lin2_max-lin2_min)*k/static_cast<double>(nx-1);
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

	void Grid::initGauLeg(double x1, double x2, std::vector<double>& Xi, std::vector<double>& Wi)
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
		std::vector<double> const& intervals,
		std::vector<double> const& sizes)
	{
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
			| std::views::transform([&](int i){ return std::make_tuple(i, _Xi[i], _Wi[i], _interval_sizes[i]); });
		for (auto const& [i, X, W, s] : tie) {
			auto min = intervals[i];
			auto max = intervals[i+1];
			
			log(LOG_DEBUG, "Grid::splitConvolution()", "Interval of size {} in range [{}, {}]", s, min, max);
			log(LOG_DEBUG, "Grid::splitConvolution()", "Abscissae: ");
			for (auto x : X)
				log(LOG_DEBUG, "Grid::splitConvolution()", "  - {}", x);
		}
		
		assert(_Xi.size() == _Wi.size() && _Xi.size() == _interval_sizes.size(), "Failed to split convolution intervals.");
		_split_interval = true;
	}

    int Grid::interpFindIdx(double x)
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

	double Grid::interpolate(ArrayGrid& yy, double x)
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
				if (std::abs(ho-hp) < 1e-15)
					log(LOG_ERROR, "Grid::interpolate()", "found a denominator equal to 0.0 ({}, {}, {})", x, xa[i], xa[i+m]);

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
		auto logx = params->logx;
		auto eplus1 = params->eplus1;
		auto& A = params->A;
		auto& E = params->E;
		auto print_count = params->print_count;

		double a = std::pow(x, 1.0-y);
		double b = std::pow(x, y);

		double interp1 = g.interpolate(A, b);
		double interp2 = g.interpolate(A, a);

		double erega = E.calcRegular(a);
		double eplusb = E.calcPlus(b);

		double res = -logx*a*erega*interp1;
		res -= logx*b*(eplusb*interp2 - eplus1*A[k])/(1.0-b);

		return res;
	}

	double Grid::convolution(ArrayGrid& A, Expression &E, uint k)
	{
		static int print_count = 0;

		double x = _points[k];
		double logx =  std::log(x);
		double eplus1 = E.plus(1.0);
		double ed1 = E.delta(1.0);
		double res = (eplus1*std::log1p(-x) + ed1) * A[k];

		auto gauleg_conv = [&](gauleg_type const& X, gauleg_type const& W, uint s) {
			for (uint i=0; i<s; i++) {
				double y = X[i];
				double w = W[i];

				double a = std::pow(x, 1.0-y);
				double b = std::pow(x, y);

				double interp1 = interpolate(A, b);
				double interp2 = interpolate(A, a);

				double erega = E.regular(a);
				double eplusb = E.plus(b);

				res -= w*logx*a*erega*interp1;
				res -= w*logx*b*(eplusb*interp2 - eplus1*A[k])/(1.0-b);
			}
		};
		
		if (!_split_interval) {
			gauleg_conv(_Xi[0], _Wi[0], _interval_sizes[0]);
		} else {
			if (!_use_gsl_routine_for_high_x) {
			    for (uint i=0; i<_Xi.size(); ++i) {
					auto const& X = _Xi[i];
					auto const& W = _Wi[i];
					auto s = _interval_sizes[i];
					gauleg_conv(X, W, s);
				}
			} else {
				for (uint i=0; i<_Xi.size()-1; ++i) {
					auto const& X = _Xi[i];
					auto const& W = _Wi[i];
					auto s = _interval_sizes[i];
					gauleg_conv(X, W, s);
				}

				GSLIntegrationParams p{
					.g = *this,
					.x = x,
					.k = k,
					.logx = logx,
					.eplus1 = eplus1,
					.A = A,
					.E = E,
					.res = ConvolutionRes{}};
				
				gsl_function f{
					.function = gsl_convolution_function,
					.params = reinterpret_cast<void*>(&p) };
				double a = 0.8, b = 1.0-1e-10;
				double epsabs = 0.0, epsrel = 1e-5;
				int key = GSL_INTEG_GAUSS41;
				std::size_t limit = gsl::DEFAULT_WORKSPACE_SIZE;
				double out{}, abserr{};
				gsl_integration_workspace* w = _workspaces[thread_index+1].get();

				int rc = gsl_integration_qags(
					&f, a, b, epsabs, epsrel,
					limit, w,
					&out, &abserr);
				p.res.out = out;
				if (rc != GSL_SUCCESS) {
					if (print_count++ > 0) 
						log(LOG_WARNING, "Grid::convolution()", "GSL integration routine failed for x = {: }", print_count, x);
					if (print_count == 0)
						log(LOG_WARNING, "Grid::convolution()", "Suppressing further GSL relating warnings.");
				}
				
				res += out;
			}
		}
		
		return res;
	}
}
