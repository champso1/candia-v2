#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"

#include <print>
#include <cmath>

namespace Candia2
{
	void Expression::fill(array_type const& grid_points, array_type const& gauss_points)
	{
		// no matter what, the plus and delta distributions
		// are evaluated at 1
		_plus_cache[1.0] = _plus_func(1.0);
		_delta_cache[1.0] = _delta_func(1.0);

		for (double x : grid_points)
		{
			for (double y : gauss_points)
			{
				double a = std::pow(x, 1.0-y);
				double b = std::pow(x, y);
			
				_reg_cache[a] = _reg_func(a);
				_plus_cache[b] = _plus_func(b);
			}
		}
	}


    double Expression::operator()(double x, uint function_part)
	{
		switch(function_part)
		{
			case REGULAR: return regular(x);
			case PLUS: return plus(x);
			case DELTA: return delta(x);
		}
		std::println("[Expression: ERROR] operator()(): Invalid function part ({}).", function_part);
		exit(EXIT_FAILURE);
	}
	
	double Expression::_reg_func(double x) const
	{
		UNUSED(x);
		return 0.0;
	}

	double Expression::_plus_func(const double x) const
	{
		UNUSED(x);
		return 0.0;
	}

	double Expression::_delta_func(const double x) const
	{
		UNUSED(x);
		return 0.0;
	}
	
	
}; // namespace Candia2
