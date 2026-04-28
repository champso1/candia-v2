/**
 *  @file Expression.hpp
 *  @brief Contains the @a Expression class and its derivations that implement functions with a regular, plus, and delta piece.
 */

#pragma once

#include <unordered_map>
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
		using cache_type = std::unordered_map<double, double>; //!< alias for cache
		using array_type = std::vector<double>; //!< alias for passed-in grid array
		using mapping_type = std::function<std::pair<double,double>(double,double)>; //!< alias for mappings
		
	protected:
		cache_type _reg_cache{}; //!< stores the values of the regular part of the expression
		cache_type _plus_cache{}; //!< stores the values of the plus part of the expression
		double _delta_cache; //!< stores the value of the delta coefficient

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
			_delta_cache = 0.0;
		}

		/**
		 *  @brief Fills the cache with the gauss-legendre points given @a grid_points and @a gauss_points
		 *  @param grid_points The array of grid points.
		 *  @param gauss_points The array of gauss-legendre abscissae
		 */
		virtual void fill(
			array_type const& grid_points, array_type const& gauss_points,
			std::span<mapping_type> const& mapping);


		/**
		 *  @defgroup piececalulators Expression Calculators
		 *  @{
		 */
		inline virtual double calcRegular([[maybe_unused]] double x) const { return 0.0; }
		inline virtual double calcPlus([[maybe_unused]] double x) const { return 0.0; }
		inline virtual double calcDelta() const { return 0.0; }
		/** @} */

		/**
		 *  @defgroup pieceretrievers Expression Retrievers
		 *  @{
		 */
		inline virtual double regular(double x) { return _reg_cache[x]; }
		inline virtual double plus(double x) { return _plus_cache[x]; }
		inline virtual double delta() { return _delta_cache; }
		/** @} */

		/** @brief calculates all nf-dependent/constant pieces of an expression */
		inline virtual void preCalc() {}
	};
};
