/**
 *  @file Expression.hpp
 *  @brief Contains the @a Expression class and its derivations that implement functions with a regular, plus, and delta piece.
 */

#pragma once

#include "Candia-v2/Common.hpp"
#include "Candia-v2/ArrayGrid.hpp"

#include <vector>

namespace Candia2
{
	enum class ExprName : uint
	{
		P0ns=0, P0qq, P0qg, P0gq, P0gg,
		P1nsm, P1nsp, P1qq, P1qg, P1gq, P1gg,
		P2nsm, P2nsp, P2nsv, P2qq, P2qg, P2gq, P2gg,
		P3nsm, P3nsp, P3nsv, P3qq, P3qg, P3gq, P3gg,
		A2ns, A2hq, A2hg, A2gq, A2gg,
		A3nsm, A3nsp, A3gq, A3gg, A3hq, A3hg, A3psqq, A3sqg, A3PSshq,
		P0ff, P0fy, P0yf, P0yy,
		Count
	};

	enum class P3ApproxType : uint
	{
		Imod1,
		Imod2,
		ImodAvg,
		Exact
	};
	
	/**
	 *  @brief Main base class that defines default behavior every "three-part" expression should have.
	 */
	class Expression
	{
	public:
		using cache_type = ArrayGrid; //!< alias for cache
		using array_type = std::vector<double>; //!< alias for passed-in grid array

	protected:
		double _plus_cache{}; //!< stores the value of the plus part of the expression
		double _delta_cache{}; //!< stores the value of the delta coefficient

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
			_plus_cache = 0.0;
			_delta_cache = 0.0;
		}

		/**
		 *  @brief Fills the cache(s) with the values of the regular part of the expression on the grid for interpolation
		 *  @param grid_points The array of grid points.
		 *  @note this function currently does nothing at the use for a regular cache has not yet manifested
		 */
		inline virtual void fill(array_type const& grid_points)
		{
		    log(LOG_WARNING, "SplittingFunction::fill()", "Not used ATM -- this method does nothing");
		}

		/**
		 *  @defgroup pieceretrievers Expression Retrievers
		 *  @{
		 */
		inline virtual double plus() { return _plus_cache; }
		inline virtual double delta() { return _delta_cache; }
		/** @} */

		/**
		 *  @defgroup piececalculates Expression Calculators
		 *  @{
		 */
		inline virtual double calcRegular([[maybe_unused]] double x) const { return 0.0; }
		inline virtual double calcPlus() const { return 0.0; }
		inline virtual double calcPlus([[maybe_unused]] double x) const { return calcPlus(); }
		inline virtual double calcDelta() const { return 0.0; }
		/** @} */
		
		/** @brief calculates all nf-dependent/constant/x-independent pieces of an expression */
		virtual void preCalc()
		{
			_plus_cache = calcPlus();
			_delta_cache = calcDelta();
		}
	};
};
