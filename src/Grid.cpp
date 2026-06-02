// Grid.cpp

#include "Candia-v2/Grid.hpp"
#include "Candia-v2/Common.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Candia2
{
    Grid::Grid(xtab_type const& xtab, GridFillerBase& grid_filler, ConvIntArgs const& convint_args)
		: _filler(grid_filler),
		  _xtab{xtab},
		  _convint_args{convint_args},
		  _Xi(convint_args.num_gauss_points, 0.0),
		  _Wi(convint_args.num_gauss_points, 0.0)
	{
		if (!std::ranges::is_sorted(xtab)) {
			log(LOG_WARNING, "Grid", "xtab array was not sorted. will sort it and continue");
			std::ranges::sort(_xtab);
		}

		grid_filler.xtab = xtab;
		grid_filler.fill(_points);
		_ntab = grid_filler.ntab;

		initGauLeg(grid_filler.min, 1.0, _Xi, _Wi);
	}

	void Grid::initGauLeg(value_type x1, value_type x2, gauleg_type& Xi, gauleg_type& Wi)
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

	constexpr auto NN = 2*INTERP_POINTS;
    int Grid::interpFindIdx(value_type x)
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

	Grid::value_type Grid::interpolate(std::span<double> yy, value_type x)
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
		uint k, value_type x, auto&& yandjaccessor,
		[[maybe_unused]] Expression& E, std::span<double> A,
		value_type plus,
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
			double fac2 = w*J * (1.0/(1.0-y))*plus*(interpa - ak);
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
		double plus = E.plus();
		double delta = E.delta();
		double res = (plus*std::log1p(-x) + delta) * A[k];

		auto mappings = _filler.get().getMappings(x);
		for (auto& mapping : mappings)
			res += mappingFunctionBase(k, x, mapping, E, A, plus, _Xi, _Wi);
		
		return res;
	}

	double Grid::convolution(std::span<double> A1, std::span<double> A2, double tau)
	{
		auto mappings = _filler.get().getMappings(tau);
		double res = 0.0;
		for (auto& mapping : mappings)
			res += mappingFunctionBase(tau, mapping, A1, A2, _Xi, _Wi);
		
		return res;
	}
}
