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
	class Grid
	{
	protected:
		std::vector<double> _grid_points; //!< the main list of grid points

		/** @name Gauss-Legendre abscissae and weights
		*/
		///@{
		std::vector<double> _Xi;
		std::vector<double> _Wi;
		///@}
		
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
		 */
		Grid(std::vector<double> const& xtab, const uint nx);

		~Grid() = default; //!< default destructor
		///@}

		/** @brief access the @a idx'th grid point (const reference)
		 *  @param idx: index of grid point
		 *  @return the @a idx'th grid point
		 */
		double const& operator[](const uint idx) const;

		/** @brief access the @a idx'th grid point (reference)
		 *  @param idx: index of grid point
		 *  @return the @a idx'th grid point
		 */
		double& operator[](const uint idx);

		/** @brief Same as the [] operator (const reference).
		 */
		inline double const& At(const uint idx) const { return operator[](idx);}

		/** @brief Same as the [] operator (reference).
		 */
		inline double& At(const uint idx) { return operator[](idx); }
		

		/** @brief Gets the size of the grid
		 *  @returns number of grid points
		 */
		inline uint Size() const { return _grid_points.size(); }


		/** @brief Returns a regular pointer to the underlying c-array
		 */
		inline double const* Ptr() const { return _grid_points.data(); }

		/** @name Setters/getters for abscissae/weights
		 */
		///@{
		inline std::vector<double> const& Abscissae() const { return _Xi; }
		inline double Abscissae(const uint idx) const { return _Xi[idx]; }
		inline std::vector<double> const& Weights() const { return _Wi; }
		inline double Weights(const uint idx) const { return _Wi[idx]; }
		///@}

		
		/** @brief Determines y(x) on the grid.
		 */
		double Interpolate(std::vector<double> const& y, const double x, bool debug=false) const;


	    /** Performs a convolution between a generic function
		 *  (represented by an array of points @a A)
		 *  and a splitting function @a P, at point k on the grid.
		 *
		 *  @param A: array of points representing the generic function
		 *  @param P: reference to @a SplittingFunction object
		 *  @param k: grid index to compute the convolution at
		 */
		double Convolution(std::vector<double> const& A,
						   std::shared_ptr<Expression> P,
						   uint k);


	private:

		/** @name Constructor helper functions
		 */
		///@{

		void InitGrid(std::vector<double> const& xtab, const uint nx); //!< fills the grid
		void InitGauLeg(); //!< init gauss-legendre abscissae/weights
		
		///@}

		
		/** returns the index pointing to the start of the range
		 *  in which to interpolate
		 */
		uint InterpFindIdx(double x) const;
	};

}


#endif // __GRID_HPP
