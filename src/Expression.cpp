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

		for (value_type x : grid_points) {
			for (value_type y : gauss_points) {
				value_type a = std::pow(x, 1.0-y);
				value_type b = std::pow(x, y);
			
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

		for (value_type x : grid_points) {
			for (value_type z : gauss_points) {
				value_type d1 = 1.0-x;
				value_type d2 = x*d1;
				value_type b = x+d1*z;
				value_type a = x/b;
			
				_reg_cache[a] = calcRegular(a);
				_plus_cache[b] = calcPlus(b);
			}
		}
	}
}; // namespace Candia2
