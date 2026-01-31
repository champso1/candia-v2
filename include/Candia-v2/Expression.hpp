/**
 *  @file Expression.hpp
 *  @brief Contains the @a Expression class and its derivations that implement functions with a regular, plus, and delta piece.
 */

#ifndef __EXPRESSION_HPP
#define __EXPRESSION_HPP

#include <map>
#include <vector>

#include "Candia-v2/Common.hpp"

namespace Candia2
{
	class Grid;

	/**
	 *  @brief Main base class that defines default behavior every "three-part" expression should have.
	 */
	class Expression
	{
	public:
		using cache_type = std::map<double, double>; //!< alias for cache
		using array_type = std::vector<double>; //!< alias for passed-in grid array
		
	protected:
		cache_type _reg_cache{}; //!< stores the values of the regular part of the expression
		cache_type _plus_cache{}; //!< stores the values of the plus part of the expression
		cache_type _delta_cache{}; //!< stores the values of the delta part of the expression

		Expression() = default; //!< default constructor
	public:
		virtual ~Expression() = default; //!< default destructor

		/** Enum to specify which part of the expression one wants */
		enum FunctionPart : uint
		{
			REGULAR,
			PLUS,
			DELTA
		};

		/**
		 *  @brief Fills the cache with the gauss-legendre points given @a grid_points and @a gauss_points.
		 *  @param grid_points The array of grid points.
		 *  @param gauss_points The array of gauss-legendre abscissae
		 */
		virtual void fill(array_type const& grid_points, array_type const& gauss_points);

		/** @brief Retrieves the regular part of the expression evaluated at x from the cache */
		inline virtual double regular(double x) { return _reg_cache[x]; }
		/** @brief Retrieves the plus part of the expression evaluated at x from the cache */
		inline virtual double plus(double x) { return _plus_cache[x]; }
		/** @brief Retrieves the delta part of the expression evaluated at x from the cache */
		inline virtual double delta(double x) { return _delta_cache[x]; }

		/**
		 *  @brief evaluates the @a function_part of the splitting function at @a x
		 */
		virtual double operator()(double x, uint function_part);
		
	protected:
		virtual double _reg_func(double x) const; //!< actually calculates the regular part of the expression at x
		virtual double _plus_func(double x) const; //!< actually calculates the plus part of the expression at x
		virtual double _delta_func(double x) const; //!< actually calculates the delta part of the expression at x
	};
	
};

#endif
