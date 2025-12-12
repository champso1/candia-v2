#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"

#include <print>

namespace Candia2
{
	double Expression::regular(double x)
	{
		if (_reg_cache.find(x) == _reg_cache.end())
		{
			double res = _reg_func(x);
			_reg_cache.emplace(x, res);
			return res;
		}
		return _reg_cache[x];
	}

	double Expression::plus(double x)
	{
		if (_plus_cache.find(x) == _plus_cache.end())
		{
			double res = _plus_func(x);
			_plus_cache.emplace(x, res);
			return res;
		}
		return _plus_cache[x];
	}

	double Expression::delta(double x)
	{
		if (_delta_cache.find(x) == _delta_cache.end())
		{
			double res = _delta_func(x);
			_delta_cache.emplace(x, res);
			return res;
		}
		return _delta_cache[x];
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
