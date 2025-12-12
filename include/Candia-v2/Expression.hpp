#ifndef __EXPRESSION_HPP
#define __EXPRESSION_HPP

#include "Candia-v2/Common.hpp"
#include <map>

namespace Candia2
{
	class Grid;
	
	class Expression
	{
	public:
		using cache_type = std::map<double, double>;
		
	protected:
		cache_type _reg_cache{}, _delta_cache{}, _plus_cache{};
		
		Expression() = default;
	public:
		virtual ~Expression() = default;

		enum FunctionPart : uint
		{
			REGULAR,
			PLUS,
			DELTA
		};

		double regular(double x);
		double plus(double x);
		double delta(double x);

		double operator()(double x, uint function_part);
		
	protected:
		virtual double _reg_func(double x) const;
		virtual double _plus_func(double x) const;
		virtual double _delta_func(double x) const;
	};
	
};

#endif
