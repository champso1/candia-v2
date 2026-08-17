// Grid.cpp

#include "Candia-v2/Grid.hpp"
#include "Candia-v2/Common.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace Candia2
{
    Grid::Grid(xtab_type const& xtab, GridFillerArgs const& gridfiller_args, ConvIntArgs const& convint_args)
		: _xtab{xtab},
		  _gridfiller_args{gridfiller_args},
		  _convint_args{convint_args},
		  _Xi(convint_args.num_gauss_points, 0.0),
		  _Wi(convint_args.num_gauss_points, 0.0)
	{
		if (!std::ranges::is_sorted(xtab)) {
			log(LOG_WARNING, "Grid", "xtab array was not sorted. will sort it and continue");
			std::ranges::sort(_xtab);
		}

		fillPoints();
		initGauLeg(gridfiller_args.min, 1.0, _Xi, _Wi);
		setupMappings();
	}

	void Grid::addXtab()
	{
		// if the number of specified points and the grid filler algorithm align just right,
		// we get an added xtab value of like 0.5 with also a 0.5000000000002 or something
		// which fucks up the interpolation algorithm
		// this just brute force checks all of the values in filled points array
		// and just removes anything that is < 1e-7 from the xtab array point
		// probably a better way to do it but who cares
		for (double x : _xtab) {
			auto val1 = std::ranges::lower_bound(_points, x);
			if (std::abs(*val1 - x) < 1e-7)
				_points.erase(val1);
			auto val2 = std::ranges::upper_bound(_points, x);
			if (std::abs(*val2 - x) < 1e-7)
				_points.erase(val2);
		}
		
		std::set<double> points_set(_points.begin(), _points.end());
		points_set.insert(_xtab.begin(), _xtab.end());
		_points = std::vector<double>(points_set.begin(), points_set.end());
		
		_ntab.clear();
		for (double x : _xtab) {
			if (auto it = std::find(_points.begin(), _points.end(), x); it != _points.end()) {
				_ntab.emplace_back(std::distance(_points.begin(), it));
				continue;
			}
		}
	}

	void Grid::fillPoints()
	{
		double log_min = std::log10(_gridfiller_args.min);
		double log_max = std::log10(_gridfiller_args.pivot1);
		uint num_log_intervals = std::round(log_max-log_min);
		double dlog = (log_max-log_min)/static_cast<double>(num_log_intervals);
		uint log_interval_size = _gridfiller_args.log_size/num_log_intervals;

		_points.clear();
		for (uint i=0; i<num_log_intervals; ++i) {
			double l0 = log_min + i*dlog;
			double l1 = l0 + dlog;
			for (uint k=0; k<log_interval_size; ++k) {
				double l = l0 + (l1-l0)*k/static_cast<double>(log_interval_size);
				_points.emplace_back(std::pow(10, l));
			}
		}

		double lin_min = _gridfiller_args.pivot1;
		double lin_max = _gridfiller_args.pivot2;
		
		for (uint k=0; k<_gridfiller_args.lin_size; ++k) {
		    double x = lin_min + (lin_max-lin_min)*k/static_cast<double>(_gridfiller_args.lin_size);
			_points.emplace_back(x);
		}

		double quad_min = _gridfiller_args.pivot2;
		double quad_max = 1.0;

		for (uint k=0; k<_gridfiller_args.quad_size; ++k) {
			double f = 1.0 - k/static_cast<double>(_gridfiller_args.quad_size);
		    double x = quad_max - (quad_max-quad_min)*f*f;
			_points.emplace_back(x);
		}
		
	    addXtab();
	}

	void Grid::initGauLeg(double x1, double x2, gauleg_type& Xi, gauleg_type& Wi)
	{
		const double eps = 3.0e-11; // relative precision

		// abscissae are symmetric:
		uint n = Xi.size(); // simpler to type
		double N = static_cast<double>(n);
		uint m = (n+1)/2;
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

	void Grid::setupMappings()
	{
		double pivot1 = _gridfiller_args.pivot1;
		double pivot2 = _gridfiller_args.pivot2;
		_mappings = {
			[=](double x, double z) {
				auto a = pivot1/x;
				return std::make_pair(x*std::pow(a, z), x*std::pow(a, z)*std::log(a));},
			[=]([[maybe_unused]] double x, double z) { return std::make_pair(pivot1+(pivot2-pivot1)*z, pivot2-pivot1); },
			[=]([[maybe_unused]] double x, double z) { return std::make_pair(1.0-(1.0-pivot2)*(1.0-z)*(1.0-z), 2.0*(1.0-pivot2)*(1.0-z)); },
			[=](double x, double z) { return std::make_pair(x+(pivot2-x)*z, pivot2-x); },
			[=]([[maybe_unused]] double x, double z) { return std::make_pair(1.0-(1.0-pivot2)*(1.0-z)*(1.0-z), 2.0*(1.0-pivot2)*(1.0-z)); },
			[=](double x, double z) { return std::make_pair(1.0-(1.0-x)*(1.0-z)*(1.0-z), 2.0*(1.0-x)*(1.0-z)); },
		};
	}

    std::span<Grid::mapping_function_type> Grid::getMappings(double x)
	{
		if (x < 0) // sanity check
			return std::span(_mappings.begin(), _mappings.end());
		else if (x > 0 && x < _gridfiller_args.pivot1)
			return std::span(_mappings.begin(), _mappings.begin()+3);
		else if (x >= _gridfiller_args.pivot1 && x < _gridfiller_args.pivot2)
			return std::span(_mappings.begin()+3, _mappings.begin()+3+2);
		else
			return std::span(_mappings.begin()+3+2, _mappings.end());
	}

	constexpr auto NN = 2*INTERP_POINTS;
    int Grid::interpFindIdx(double x)
	{
		static const int n = static_cast<int>(size());
		static const int max_k = static_cast<int>(n-NN);

		auto it = std::upper_bound(
			_points.begin(),
			_points.end(),
			x);
		int k = static_cast<int>(it - _points.begin()) - INTERP_POINTS;

		return std::clamp(k, 0, max_k);
	}

    double Grid::interpolate(std::span<double> yy, double x)
	{
		constexpr int n = 2 * INTERP_POINTS;
    
		int k = interpFindIdx(x);
		double const* xa = &(_points[k]);
		double const* ya = &(yy[k]);

		double c[n];
		double d[n];
    
		int ns = 0;
		double dif = std::abs(x - xa[0]);
		c[0] = ya[0];
		d[0] = ya[0];

		for (int i=1; i<n; i++) {
			double dift = std::abs(x - xa[i]);
			if (dift < dif) {
				ns = i;
				dif = dift;
			}
			c[i] = ya[i];
			d[i] = ya[i];
		}

		double y = ya[ns--];

		for (int m=1; m<n; m++) {
			for (int i=0; i<(n-m); i++) {
				double ho = xa[i] - x;
				double hp = xa[i+m] - x;
				double w = c[i+1] - d[i];

				double grid_diff = xa[i] - xa[i+m]; 

				if (grid_diff > -1e-15) {
					log(LOG_ERROR, "Grid::interpolate()", "found a denominator equal to 0.0 ({}, {}, {})", x, xa[i], xa[i+m]);
				}

				double den = w/grid_diff;
				d[i] = hp*den;
				c[i] = ho*den;
			}

			y += (2*ns < (n-1-m) ? c[ns+1] : d[ns--]);
		}

		return y;
	}

	double Grid::mappingFunctionBase(
		uint k, double x, auto&& yandjaccessor,
		Expression& E, std::span<double> A,
		double plus,
		gauleg_type const& X, gauleg_type const& W)
	{
		double ak = A[k];
		double out = 0.0;
		uint s = X.size();
		for (uint i=0; i<s; i++) {
			double z = X[i];
			double w = W[i];
			auto [y, J] = yandjaccessor(x, z);
			double a = x/y;
			double eregy = E.calcRegular(y);
			double interpa = interpolate(A, a);

			double fac1 = w*J * eregy*interpa;
			double fac2 = w*J * (1.0/(1.0-y))*(E.calcPlus(y)*interpa - plus*ak);
			out += fac1;
			out += fac2;
		}
		return out;
	}
	double Grid::mappingFunctionBase(
	    double tau, auto&& yandjaccessor,
		std::span<double> A1, std::span<double> A2,
		gauleg_type const& X, gauleg_type const& W)
	{
		double out = 0.0;
		for (uint i=0; i<X.size(); i++) {
			double z = X[i];
			double w = W[i];
			auto [x, J] = yandjaccessor(tau, z);
			double a = tau/x;
			double interp_a1 = interpolate(A1, x);
			double interp_a2 = interpolate(A2, a);

			out += w*J * (1/x) * interp_a1*interp_a2;
		}
		return out;
	}

	double Grid::convolution(std::span<double> A, Expression &E, uint k)
	{
		double x = _points[k];
		double plus1 = E.calcPlus(1.0);
		double delta = E.delta();
		double res = (plus1*std::log1p(-x) + delta) * A[k];

		auto mappings = getMappings(x);
		for (auto& mapping : mappings)
			res += mappingFunctionBase(k, x, mapping, E, A, plus1, _Xi, _Wi);
		
		return res;
	}

	double Grid::convolution(std::span<double> A1, std::span<double> A2, double tau)
	{
		auto mappings = getMappings(tau);
		double res = 0.0;
		for (auto& mapping : mappings)
			res += mappingFunctionBase(tau, mapping, A1, A2, _Xi, _Wi);
		
		return res;
	}
}
