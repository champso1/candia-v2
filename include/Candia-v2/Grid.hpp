/** @file
 *
 *  Contains the information for the grid points on which the evolution is run.
 *  Also contains functionality for convolution with the splitting functions
 *  along with basic polynomial interpolation.
 */

#ifndef __GRID_HPP
#define __GRID_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"

#include <memory>
#include <vector>

namespace Candia2
{
	
	/** @brief Class to handle the logarithmic grid
	 */
	class Grid final
	{
	protected:
		using grid_type = std::vector<double>;
		using gauleg_type = std::vector<double>;
		using ntab_type = std::vector<int>;
		
		grid_type _points{};
		
		/** @name Gauss-Legendre abscissae and weights
		*/
		///@{
		gauleg_type _Xi{};
		gauleg_type _Wi{};
		///@}

		/** @name additional G-L abscissae and weights
		*/
		///@{
		gauleg_type _Xi_low{}, _Xi_high{};
		gauleg_type _Wi_low{}, _Wi_high{};
		///@}	

		ntab_type _ntab; //!< stores indices for the tabulated grid points
	public:

		/** @name Constructors/destructors
		 */
		///@{
		/** @brief Deleted default constructor
		 *
		 *  Must provide some argument to initialize the grid.
		 */
		Grid() = delete;

		/** @param xtab: list of tabulated x-values.
		 *  @param nx: number of desired grid points
		 *  @note Only fills grid points, doesn't evaluate any function
		 */
		Grid(std::vector<double> const& xtab, const uint nx);

		/** @brief default destructor */
		~Grid() = default;
		///@}

		/** @name Accessors
		 */
		///@{
		inline grid_type const& Points() const { return _points; }
		inline double At(const uint idx) const { return _points.at(idx); };
		inline double operator[](const uint idx) const { return _points.at(idx); }
		///@}

		
		/** @name Setters/getters for abscissae/weights
		 */
		///@{
		inline gauleg_type const& Abscissae() const { return _Xi; }
		inline double Abscissae(uint idx) const { return _Xi[idx]; }
		inline gauleg_type const& Weights() const { return _Wi; }
		inline double Weights(uint idx) const { return _Wi[idx]; }
		///@}

		inline uint Size() const { return _points.size(); }

		/** @name Setters/getters for @a ntab array
		 */
		///@{
		inline ntab_type const& Ntab() const { return _ntab; }
		inline int const& Ntab(uint idx) const { return _ntab.at(idx); }
		inline ntab_type& Ntab() { return _ntab; }
		inline int& Ntab(uint idx) { return _ntab.at(idx); }
		///@}

		/** @name c++ std range-based functions
		 */
		///@{
		inline grid_type::const_iterator begin() const { return _points.begin(); }
		inline grid_type::iterator begin() { return _points.begin(); }
		inline grid_type::const_iterator end() const { return _points.end(); }
		inline grid_type::iterator end() { return _points.end(); }
		///@}

		
		/** @brief Determines y(x) on the grid
		 *
		 *  @param y: vector of points of the function's values on the grid
		 *  @param x: the point to interpolate to
		 *
		 *  @note This is only for functions stored within a generic vector
		 *  not for functions stored on a function grid. 
		 */
		double Interpolate(grid_type const& y, const double x) const;


	    /** Performs a convolution between a generic function
		 *  (represented by an array of points @a A)
		 *  and a splitting function @a P, at point k on the grid.
		 *
		 *  @param A: array of points representing the generic function
		 *  @param P: reference to @a SplittingFunction object
		 *  @param k: grid index to compute the convolution at
		 *  @param highorder_conv: whether to split convolution into high (>0.9) and low regions
		 *
		 *  @return The value of the convolution
		 *
		 *  @note This is only for functions evaluated on a generic vector,
		 *  if a function is evaluated on a function grid, use that function's Interpolate
		 */
		double Convolution(grid_type const& A,
						   std::shared_ptr<Expression> P,
						   uint k, bool highorder_conv=false);

	private:

		/** @name Constructor helper functions
		 */
		///@{
		void InitGrid(grid_type const& xtab, const uint nx); //!< fills the grid
		void InitGrid2(grid_type const& xtab, const uint nx); //!< fills the grid
		void InitGrid3(grid_type const& xtab, const uint nx); //!< fills the grid 
		void InitGauLeg(double x1, double x2, std::vector<double> & Xi, std::vector<double> & Wi); //!< init gauss-legendre abscissae/weights
		///@}


	public:
		/** returns the index pointing to the start of the range
		 *  in which to interpolate for the value @a x
		 */
		uint InterpFindIdx(double x) const;
	};
}

#endif // __GRID_HPP
