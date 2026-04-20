// Expression.cpp

#include "Candia-v2/Expression.hpp"

namespace Candia2
{
    void Expression::fill(
		array_type const& grid_points, array_type const& gauss_points,
		std::span<mapping_type> const& mappings)
	{
		// no matter what, the plus and delta distributions
		// are evaluated at 1
		_plus_cache[1.0] = calcPlus(1.0);
		_delta_cache = calcDelta();

		for (double x : grid_points) {
			for (double z : gauss_points) {
				for (auto& mapping : mappings) {
					auto [y,_] = mapping(x, z);
			
					_reg_cache[y] = calcRegular(y);
					_plus_cache[y] = calcPlus(y);
				}
			}
		}
	}
} // namespace Candia2
