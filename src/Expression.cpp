#include "Candia-v2/Expression.hpp"

#include <cmath>

namespace Candia2
{
	void Expression::fill(array_type const& grid_points, array_type const& gauss_points)
	{
		// no matter what, the plus and delta distributions
		// are evaluated at 1
		_plus_cache[1.0] = calcPlus(1.0);
		_delta_cache[1.0] = calcDelta(1.0);

		for (double x : grid_points) {
			for (double y : gauss_points) {
				double a = std::pow(x, 1.0-y);
				double b = std::pow(x, y);
			
				_reg_cache[a] = calcRegular(a);
				_plus_cache[b] = calcPlus(b);
			}
		}
	}

	void Expression::fill2(array_type const& grid_points, array_type const& gauss_points)
	{
		// no matter what, the plus and delta distributions
		// are evaluated at 1
		_plus_cache[1.0] = calcPlus(1.0);
		_delta_cache[1.0] = calcDelta(1.0);

		for (double x : grid_points) {
			for (double z : gauss_points) {
				double d1 = 1.0-x;
				double d2 = x*d1;
				double b = x+d1*z;
				double a = x/b;
			
				_reg_cache[a] = calcRegular(a);
				_plus_cache[b] = calcPlus(b);
			}
		}
	}

    double Expression::operator()(double x, uint function_part)
	{
		switch(function_part) {
			case REGULAR: return regular(x);
			case PLUS: return plus(x);
			case DELTA: return delta(x);
		}
		log(LOG_ERROR, "Expression::operator()()", "Invalid function part ({}).", function_part);
		return 0.0; // unreachable
	}
	
}; // namespace Candia2
