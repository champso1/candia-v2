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

		/** @brief clears all the cache values */
		inline virtual void clear()
		{
			_reg_cache.clear();
			_plus_cache.clear();
			_delta_cache.clear();
		}

		/**
		 *  @brief Fills the cache with the gauss-legendre points given @a grid_points and @a gauss_points in the small-x mapping
		 *  @param grid_points The array of grid points.
		 *  @param gauss_points The array of gauss-legendre abscissae
		 */
		virtual void fill(array_type const& grid_points, array_type const& gauss_points);

		/**
		 *  @brief Fills the cache with the gauss-legendre points given @a grid_points and @a gauss_points in the large-x mapping
		 *  @param grid_points The array of grid points.
		 *  @param gauss_points The array of gauss-legendre abscissae
		 */
		virtual void fill2(array_type const& grid_points, array_type const& gauss_points);

		/**
		 *  @brief Fills the cache with the gauss-legendre points given @a grid_points and @a gauss_points.
		 *
		 *  This method takes in multiple sets of gauss_points, which will be the case
		 *  if the user splits the grid into convolution intervals.
		 *
		 *  @param grid_points The array of grid points.
		 *  @param gauss_points The array of several sets of gauss-legendre abscissae
		 */
		virtual inline void fill(array_type const& grid_points, std::vector<array_type> const& gauss_points)
		{
			clear();
			for (auto it = gauss_points.begin(); it!=gauss_points.end(); ++it)
				fill(grid_points, *it);
		}

		/** @brief actually calculates the regular distribution */
		inline virtual double calcRegular(double x) const { return 0.0; }
		/** @brief actually calculates the plus distribution */
		inline virtual double calcPlus(double x) const { return 0.0; }
		/** @brief actually calculates the delta distribution */
		inline virtual double calcDelta(double x) const { return 0.0; }
		

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
	};
	
};

#endif
