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
	/**
	 *  @brief Main base class that defines default behavior every "three-part" expression should have.
	 */
	class Expression
	{
	public:
		using value_type = double;
		using cache_type = std::map<value_type, value_type>; //!< alias for cache
		using array_type = std::vector<value_type>; //!< alias for passed-in grid array
		using mapping_type = std::vector<std::function<std::pair<double,double>(double,double)>>; //!< alias for mappings
		
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
		 *  @brief Fills the cache with the gauss-legendre points given @a grid_points and @a gauss_points
		 *  @param grid_points The array of grid points.
		 *  @param gauss_points The array of gauss-legendre abscissae
		 */
		virtual void fill(
			array_type const& grid_points, array_type const& gauss_points,
			mapping_type::value_type const& mapping);

		/**
		 *  @brief Fills the cache with the gauss-legendre points given @a grid_points and @a gauss_points.
		 *
		 *  This method takes in multiple sets of gauss_points, which will be the case
		 *  if the user splits the grid into convolution intervals.
		 *
		 *  @param grid_points The array of grid points.
		 *  @param gauss_points The array of several sets of gauss-legendre abscissae
		 */
		virtual inline void fill(
			array_type const& grid_points, std::vector<array_type> const& gauss_points,
			mapping_type const& mappings)
		{
			clear();
			for (auto gauleg_it = gauss_points.begin(); gauleg_it!=gauss_points.end(); ++gauleg_it) {
				for (auto mapping_it = mappings.begin(); mapping_it != mappings.end(); ++mapping_it) {
					fill(grid_points, *gauleg_it, *mapping_it);
				}
			}
		}

		/** @brief actually calculates the regular distribution */
		inline virtual value_type calcRegular(value_type x) const { return 0.0; }
		/** @brief actually calculates the plus distribution */
		inline virtual value_type calcPlus(value_type x) const { return 0.0; }
		/** @brief actually calculates the delta distribution */
		inline virtual value_type calcDelta(value_type x) const { return 0.0; }
		

		/** @brief Retrieves the regular part of the expression evaluated at x from the cache */
		inline virtual value_type regular(value_type x) { return _reg_cache[x]; }
		/** @brief Retrieves the plus part of the expression evaluated at x from the cache */
		inline virtual value_type plus(value_type x) { return _plus_cache[x]; }
		/** @brief Retrieves the delta part of the expression evaluated at x from the cache */
		inline virtual value_type delta(value_type x) { return _delta_cache[x]; }
	};
};

#endif
