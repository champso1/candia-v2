/**
 *  @file Grid.hpp
 *  @brief Contains the @a Grid class which contains the interpolation/convolution grids and associated routines.
 */

#pragma once

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"

namespace Candia2
{
	/** @brief Struct for grouping the paramaters related to convolution/integration */
	struct GridFillerArgs final
	{
		double min{1.0e-5};
		uint log_size{100};
		uint lin_size{50};
		uint quad_size{25};
		double pivot1{0.1}, pivot2{0.9};
	};

	/** @brief Struct for grouping the paramaters related to convolution/integration */
	struct ConvIntArgs final
	{
		uint num_gauss_points{50};
		uint num_interp_points{4};
	};
	/**
	 *  @brief Class that contains the interpolation/convolution grid and methods
	 */
	class Grid final
	{
    public:
		using grid_type = std::vector<double>; //!< alias for the underlying grid type
		using gauleg_type = std::vector<double>; //!< alias for the type of the array of gauss-legendre weights/abscissae
		using ntab_type = std::vector<int>; //!< alias for the type of the calulated ntab array
		using xtab_type = std::vector<double>;
	private:
		grid_type _points; //!< grid points
		ntab_type _ntab{};     //!< stored indices for the tabulated grid points
		xtab_type _xtab{};     //!< stored values of the tabulated grid points
		GridFillerArgs _gridfiller_args; //!< contains options related to filling the grid points
		ConvIntArgs _convint_args; //!< contains misc convolution/interpolation options/args
		gauleg_type _Xi{}; //!< list of split-up gauleg abscissae per interval
		gauleg_type _Wi{}; //!< list of split-up gauleg weights per interval

		// mapping info
	private:
		using mapping_function_type = std::function<std::pair<double,double>(double,double)>;
		std::vector<mapping_function_type> _mappings;
		void setupMappings();
		std::span<mapping_function_type> getMappings(double x);
		
	public:
		/** @brief provides a facility like python's enumerate() to return an index and value at once */
		struct EnumerateIterator final
		{
			grid_type const& data;

			struct Iterator final
			{
				grid_type const& v;
				grid_type::size_type idx;

				inline auto operator*() const {
					return std::pair<uint, double>{idx, v[idx]};
				}

				inline Iterator& operator++() { ++idx; return *this; }
				inline bool operator!=(Iterator const& other) { return idx != other.idx; }
			};

			inline auto begin() { return Iterator{data, 0}; }
			inline auto end() { return Iterator{data, data.size()}; }
		};

		inline auto enumerate() {
			return EnumerateIterator{_points};
		}
		
	public:
		Grid() = delete; //!< must provide information to fill the grid
		/**
		 *  fills the grid, sets up GSL objects, and other initialization functions
		 *  @param xtab Array of tabulated grid points to ensure the grid contains for easy retrieval
		 *  @param grid_filler an object that will fill the grid in a particular way
		 *  @param gauleg_args arguments to setup/initialize how gauleg integration behaves
		 */
		Grid(xtab_type const& xtab, GridFillerArgs const& grid_filler={}, ConvIntArgs const& gauleg_args={});
		Grid(Grid const& other) = default;
		Grid(Grid&& other) = default;
		~Grid() = default;

		inline uint size() const { return _points.size(); }

		inline xtab_type& xtab() { return _xtab; }
		inline xtab_type const& xtab() const { return _xtab; }
		inline ntab_type const& ntab() const { return _ntab; }
		inline ntab_type& ntab() { return _ntab; }

		inline grid_type const& points() const { return _points; }
		inline double operator[](uint idx) const { return _points[idx]; }

		inline gauleg_type const& abscissae() const { return _Xi; }
		inline gauleg_type const& weights() const { return _Wi; }

		inline auto begin() const { return _points.begin(); }
		inline auto begin() { return _points.begin(); }
		inline auto end() const { return _points.end(); }
		inline auto end() { return _points.end(); }

		/**
		 *  @brief Handles a convolution between a splitting function or OME with simple mappings for y -> z
		 *  @param k grid index
		 *  @param x x-value at the grid index
		 *  @param yandjaccessor a @a YandJAccessor to retrieve y and the jacobian given x and z (the mapped value, a gauleg abscissa)
		 *  @param E splitting function / operator matrix elements
		 *  @param A array to convolute
		 *  @param eplus1 the constant value of the plus component of the expression evaluated at x=1
		 *  @param X the list of gauleg abscissae
		 *  @param W the list of gauleg weights
		 */
		double mappingFunctionBase(
			uint k, double x, auto&& yandjaccessor,
			Expression& E, std::span<double> A,
			double eplus1,
			gauleg_type const& X, gauleg_type const& W);

		/**
		 *  @brief Handles a convolution betweewn two arrays with simple mappings for y -> z
		 *  @param tau the other factor in the convolution (since not necessarily a grid point)
		 *  @param yandjaccessor a @a YandJAccessor to retrieve y and the jacobian given x and z (the mapped value, a gauleg abscissa)
		 *  @param A1 array1 to convolute
		 *  @param A2 array2 to convolute
		 *  @param X the list of gauleg abscissae
		 *  @param W the list of gauleg weights
		 */
		double mappingFunctionBase(
		    double tau, auto&& yandjaccessor,
		    std::span<double> A1, std::span<double> A2,
			gauleg_type const& X, gauleg_type const& W);

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
		double interpolate(std::span<double> y, double x);
		/**
		 *  @brief Performs a convolution between an array @a A and an expression @a E
		 *  @param A the array
		 *  @param E the expression
		 *  @param k the grid index to perform the convolution at
		 */
		double convolution(std::span<double> A, Expression& E, uint k);

		/**
		 *  @brief Performs a convolution between two arrays
		 *  @param A1 array1
		 *  @param A2 array2
		 *  @param q energy (squared)
		 *  @param tau the other factor used in the convolution
		 */
		double convolution(std::span<double> A1, std::span<double> A2, double tau);
	private:
		/** @brief adds the xtab array to the set of points, also filling in ntab */
		void addXtab();
		/** @brief fills the grid points in  */
		void fillPoints();
	    /** @brief Initializes the set of gauss-legendre weights and abscissae */
		void initGauLeg(double x1, double x2, gauleg_type& Xi, gauleg_type& Wi);
	};
}
