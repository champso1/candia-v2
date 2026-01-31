/**
 *  @file Grid.hpp
 *  @brief Contains the @a Grid class which contains the interpolation/convolution grids and associated routines.
 */

#ifndef __GRID_HPP
#define __GRID_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"

#include <vector>

namespace Candia2
{
	class ArrayGrid;
	/**
	 *  @brief Class that contains the interpolation/convolution grid and the methods to perform the interpolation and convolution.
	 */
	class Grid final
	{
	public:
		using grid_type = std::vector<double>; //!< alias for the underlying grid type
		using gauleg_type = std::vector<double>; //!< alias for the type of the array of gauss-legendre weights/abscissae
		using ntab_type = std::vector<int>; //!< alias for the type of the calulated ntab array

	private:
		grid_type _points{}; //!< grid points
		ntab_type _ntab;     //!< stored indices for the tabulated grid points
		grid_type _xtab;     //!< stored values of the tabulated grid points

		uint _gauss_points;  //!< the number of gauss-legendre points
		gauleg_type _Xi{};   //!< the array of gauss-legendre abscissae
		gauleg_type _Wi{};   //!< the array of gauss-legendre weight
	public:
		Grid() = delete; //!< default constructor deleted; must provide information to fill the grid
		/**
		 *  Fills the grid with @a nx grid points according to @a grid_fill_type and sets up @a gauss_points gauss-legendre points
		 *  @param xtab Array of tabulated grid points to ensure the grid contains for easy retrieval
		 *  @param nx number of grid points
		 *  @param gauss_points number of gauss_legendre points
		 *  @param grid_fill_type Specifies how the grid is laid out
		 */
		Grid(grid_type const& xtab, uint nx, uint gauss_points, int grid_fill_type=1);
		~Grid() = default; //!< default destructor

		/** Getter for xtab array */
		inline grid_type& xtab() { return _xtab; }
		/** Getter for xtab array */
		inline grid_type const& xtab() const { return _xtab; }

		/** Getter for point on the grid */
		inline grid_type const& points() const { return _points; }
		/** Getter for point on the grid */
		inline double at(uint idx) const { return _points.at(idx); };
		/** Getter for point on the grid */
		inline double operator[](uint idx) const { return _points[idx]; }

		/** Getter for the gauss-legendre abscissae */
		inline gauleg_type const& abscissae() const { return _Xi; }
		/** Getter for the gauss-legendre abscissae */
		inline double abscissae(uint idx) const { return _Xi[idx]; }
		/** Getter for the gauss-legendre weights */
		inline gauleg_type const& weights() const { return _Wi; }
		/** Getter for the gauss-legendre weights */
		inline double weights(uint idx) const { return _Wi[idx]; }

		/** Getter for the grid size */
		inline uint size() const { return _points.size(); }

		/** Getter for the ntab array */
		inline ntab_type const& ntab() const { return _ntab; }
		/** Getter for the ntab array */
		inline int const& ntab(uint idx) const { return _ntab.at(idx); }
		/** Getter for the ntab array */
		inline ntab_type& ntab() { return _ntab; }
		/** Getter for the ntab array */
		inline int& ntab(uint idx) { return _ntab.at(idx); }

		/** Const iterator to the beginning of the underlying array */
		inline grid_type::const_iterator begin() const { return _points.begin(); }
		/** Iterator to the beginning of the underlying array */
		inline grid_type::iterator begin() { return _points.begin(); }
		/** Const iterator to the end of the underlying array */
		inline grid_type::const_iterator end() const { return _points.end(); }
		/** Iterator to the end of the underlying array */
		inline grid_type::iterator end() { return _points.end(); }

		/**
		 *  @brief Uses a binary search to find the grid point closest to the given value of x
		 *  @param x value to search for
		 */
		int interpFindIdx(double x);
		/**
		 *  @brief Interpolates array @a y at @a x on the grid.
		 *  @param y the array grid to interpolate
		 *  @param x the value of x to interpolate at
		 */
		double interpolate(ArrayGrid& y, double x);
		/**
		 *  @brief Performs a convolution between an array @a A and an expression @a E
		 *  @param A the array
		 *  @param E the expression
		 *  @param k the grid index to perform the convolution at
		 */
		double convolution(ArrayGrid& A, Expression &E, uint k);
		
	private:
		/** Fills the grid according to the original candia-v2 method (log-spaced) */
		void initGrid(grid_type const& xtab, uint nx);
		/** Fills the grid with log-spacing and a calculated set of additional points from \f$0.1<x<1.0\f$ */
		void initGrid2(grid_type const& xtab, uint nx);
		/** Fills the grid with log-spacing and a pre-defined additional set of points from \f$0.1<x<1.0\f$ */
		void initGrid3(grid_type const& xtab, uint nx);
		/** Initializes the set of gauss-legendre weights and abscissae */
		void initGauLeg(double x1, double x2, std::vector<double> & Xi, std::vector<double> & Wi);
	};
}

#endif // __GRID_HPP
