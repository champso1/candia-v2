#ifndef __EXPRESSION_HPP
#define __EXPRESSION_HPP

#include "Candia-v2/Common.hpp"
#include <map>
#include <vector>

namespace Candia2
{
	class Grid;
	
	class Expression
	{
	public:
		using cache_type = std::map<double, double>;
		using array_type = std::vector<double>;
		
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

		virtual void fill(array_type const& grid_points, array_type const& gauss_points);
		// virtual void add_point(double x, uint function_part);

		inline virtual double regular(double x) { return _reg_cache[x]; }
		inline virtual double plus(double x) { return _plus_cache[x]; }
		inline virtual double delta(double x) { return _delta_cache[x]; }

		virtual double operator()(double x, uint function_part);
		
	protected:
		virtual double _reg_func(double x) const;
		virtual double _plus_func(double x) const;
		virtual double _delta_func(double x) const;
	};
	
};

#endif
