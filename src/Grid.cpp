// Grid.cpp

#include "Candia-v2/Grid.hpp"
#include "Candia-v2/Common.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Candia2
{
    Grid::Grid(grid_type const& xtab, GridFillerBase& grid_filler, ConvIntArgs const& convint_args)
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

	static thread_local int prev_interp_idx = 0;
    int Grid::interpFindIdx(value_type x)
	{
		const int n = static_cast<int>(size());
		const int max_k = static_cast<int>(n-2*INTERP_POINTS);

		auto it = std::upper_bound(_points.begin() + prev_interp_idx, _points.end(), x);
		
		int k = static_cast<int>(it - _points.begin()) - INTERP_POINTS;
		
		if (k < 0)
			k = 0;
		else if (k > max_k)
			k = max_k;

		prev_interp_idx = k;

		return k;
	}

	Grid::value_type Grid::interpolate(std::span<double> yy, value_type x)
	{
		constexpr int n = 2*INTERP_POINTS;
		int ns=0;
		double y, den, dif, dift, ho, hp, w;

		int k = interpFindIdx(x);
		
		double const* xa = &(_points.data()[k]);
		double const* ya = &(yy.data()[k]);

		double c[n]{};
		double d[n]{};
		
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

	Grid::value_type Grid::mappingFunctionBase(
		uint k, value_type x, auto&& yandjaccessor,
		Expression& E, std::span<double> A,
		value_type eplus1,
		gauleg_type const& X, gauleg_type const& W)
	{
		auto reg = [&](double y) {
			if (_use_cached_expressions)
				return E.regular(y);
			else
				return E.calcRegular(y);
		};
		auto plus = [&](double y) {
			if (_use_cached_expressions)
				return E.plus(y);
			else
				return E.calcPlus(y);
		};
		
		double ak = A[k];
		double out = 0.0;
		uint s = X.size();
		prev_interp_idx = 0;
		for (uint i=0; i<s; i++) {
			double z = X[i];
			double w = W[i];
			auto [y, J] = yandjaccessor(x, z);
			double a = x/y;
			double eregy = reg(y);
			double interpa = interpolate(A, a);
			double eplusy = plus(y);
			
			out += w*J * eregy*interpa;
			out += w*J * (1.0/(1.0-y))*(eplusy*interpa - eplus1*ak);
		}
		return out;
	}

	Grid::value_type Grid::convolution(std::span<double> A, Expression &E, uint k)
	{
		double x = _points[k];
		double eplus1 = _use_cached_expressions ? E.plus(1.0) : E.calcPlus(1.0);
		double ed1 = _use_cached_expressions ? E.delta() : E.calcDelta();
		double res = (eplus1*std::log1p(-x) + ed1) * A[k];

		auto mappings = _filler.get().getMappings(x);
		for (auto& mapping : mappings)
			res += mappingFunctionBase(k, x, mapping, E, A, eplus1, _Xi, _Wi);
		
		return res;
	}
}
