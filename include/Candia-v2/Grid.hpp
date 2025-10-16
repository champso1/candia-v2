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

#include <iterator>
#include <memory>
#include <tuple>
#include <vector>

namespace Candia2
{
	
	/** @brief Class to handle the logarithmic grid
	 */
	class Grid final
	{
	protected:
		std::vector<double> _points{};
		
		/** @name Gauss-Legendre abscissae and weights
		*/
		///@{
		std::vector<double> _Xi{};
		std::vector<double> _Wi{};
		///@}

		std::vector<int> _ntab; //!< stores indices for the tabulated grid points
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
		inline std::vector<double> const& Points() const { return _points; }
		inline double At(const uint idx) const { return _points.at(idx); };
		inline double operator[](const uint idx) const { return _points.at(idx); }
		///@}

		
		/** @name Setters/getters for abscissae/weights
		 */
		///@{
		inline std::vector<double> const& Abscissae() const { return _Xi; }
		inline double Abscissae(const uint idx) const { return _Xi[idx]; }
		inline std::vector<double> const& Weights() const { return _Wi; }
		inline double Weights(const uint idx) const { return _Wi[idx]; }
		///@}

		inline uint Size() const { return _points.size(); }

		/** @name Setters/getters for @a ntab array
		 */
		///@{
		inline std::vector<int> const& Ntab() const { return _ntab; }
		inline int const& Ntab(const uint idx) const { return _ntab.at(idx); }
		inline std::vector<int>& Ntab() { return _ntab; }
		inline int& Ntab(const uint idx) { return _ntab.at(idx); }
		///@}

		
		/** @brief Determines y(x) on the grid
		 *
		 *  @param y: vector of points of the function's values on the grid
		 *  @param x: the point to interpolate to
		 *
		 *  @note This is only for functions stored within a generic vector
		 *  not for functions stored on a function grid. 
		 */
		double Interpolate(std::vector<double> const& y, const double x) const;


	    /** Performs a convolution between a generic function
		 *  (represented by an array of points @a A)
		 *  and a splitting function @a P, at point k on the grid.
		 *
		 *  @param A: array of points representing the generic function
		 *  @param P: reference to @a SplittingFunction object
		 *  @param k: grid index to compute the convolution at
		 *
		 *  @return The value of the convolution
		 *
		 *  @note This is only for functions evaluated on a generic vector,
		 *  if a function is evaluated on a function grid, use that function's Interpolate
		 */
		double Convolution(std::vector<double> const& A,
						   std::shared_ptr<Expression> P,
						   uint k);

	private:

		/** @name Constructor helper functions
		 */
		///@{
		void InitGrid(std::vector<double> const& xtab, const uint nx); //!< fills the grid
		void InitGrid2(std::vector<double> const& xtab, const uint nx); //!< fills the grid// 
		void InitGauLeg(); //!< init gauss-legendre abscissae/weights
		///@}


	public:
		/** returns the index pointing to the start of the range
		 *  in which to interpolate for the value @a x
		 */
		uint InterpFindIdx(double x) const;
	};


	class FunctionGrid final
	{
	private:
	    Grid const& _grid;
		std::array<std::vector<double>, 3> _y{};
		
	public:
		FunctionGrid(Grid const& grid)
			:  _grid{grid},
			   _y([](const uint size) -> decltype(_y) {
				   std::array<std::vector<double>, 3> v{};
				   for (std::vector<double>& _v : v)
					   _v.resize(size);
				   return v;
			   }(grid.Size()))
		{}

		enum FunctionPart : uint
		{
			REGULAR = 0,
			PLUS,
			DELTA
		};
		
		/** @brief Determines the value of this function on the underlying grid
		 *
		 *  @param part: which part of the expression to interpolate on
		 *  @param x: the point to interpolate to
		 */
		double Interpolate(FunctionPart part, const double x) const;

		/** Performs a convolution between this function
		 *  and a splitting function @a P, at point k on the grid.
		 *
		 *  @param A: vector of points to convolute this function with
		 *  @param k: grid index to compute the convolution at
		 */
		double Convolution(std::vector<double> A, uint k) const;

	private:
		inline double X(const uint idx) const { return _grid.At(idx); }
		inline std::vector<double> const& Y(const FunctionPart part) const
		{
			switch (part)
			{
				case REGULAR: return std::get<0>(_y);
				case PLUS: return std::get<1>(_y);
				case DELTA: return std::get<2>(_y);
			}
		}
		inline double Y(const FunctionPart part, const uint idx) const
		{
			return Y(part).at(idx);
		}
	};

}

#endif // __GRID_HPP
