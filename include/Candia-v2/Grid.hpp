/**
 *  @file Grid.hpp
 *  @brief Contains the @a Grid class which contains the interpolation/convolution grids and associated routines.
 */

#ifndef __GRID_HPP
#define __GRID_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"

#include <gsl/gsl_errno.h>
#include <vector>
#include <memory>

#include <gsl/gsl_integration.h>


namespace Candia2
{
	namespace gsl
	{
		inline auto workspace_deleter = [](gsl_integration_workspace* w){ gsl_integration_workspace_free(w); };
		using workspace_deleter_type = decltype(workspace_deleter);
		using workspace_type = std::unique_ptr<gsl_integration_workspace, workspace_deleter_type>;

		inline auto make_workspace = [](uint size){ return workspace_type(gsl_integration_workspace_alloc(size), workspace_deleter); };
		static constexpr uint DEFAULT_WORKSPACE_SIZE = 1000;
		inline auto make_default_workspace = [](){return workspace_type(gsl_integration_workspace_alloc(DEFAULT_WORKSPACE_SIZE), workspace_deleter); };

		static inline int error_print_count = 20;

		extern "C" {
			static inline void error_handler(
				const char * reason, const char * file,
				int line, int gsl_errno)
			{
				if (error_print_count > 0) {
					log(LOG_ERROR_NOQUIT, "GSL", "({}:{}) {}", file, line, reason);
					error_print_count--;
				}
				if (error_print_count == 0)
					log(LOG_WARNING, "GSL", "Reached more than 20 GSL failures, suppressing additional ones.");
			}
		}

		static bool error_handler_set = (([](){ gsl_set_error_handler(error_handler); })(), true);
	}
	
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

		enum GridFillType : uint
		{
			LOG = 0, //!< simple logarithmic intervals
			LOG_LIN, //!< logarithmic intervals until 0.1, then linear until 1.0
		};

		struct ConvolutionRes final
		{
			double out{};
			std::vector<double> y{}, w{}, a{}, b{}, interp1{}, interp2{}, erega{}, eplusb{};
		};
		struct GSLIntegrationParams final
		{
			Grid& g;

			double x;
			uint k;
			double logx;
			double eplus1;

			ArrayGrid& A;
			Expression& E;

			int print_count;
			ConvolutionRes res;
		};
	private:
		grid_type _points{}; //!< grid points
		ntab_type _ntab;     //!< stored indices for the tabulated grid points
		grid_type _xtab;     //!< stored values of the tabulated grid points

		bool _split_interval{false}; //!< flag that declares if the user has split the convolution into intervals
		std::vector<gauleg_type> _Xi{}; //!< list of split-up gauleg abscissae per interval
		std::vector<gauleg_type> _Wi{}; //!< list of split-up gauleg weights per interval
		std::vector<uint> _interval_sizes{}; //!< number of points per interval

		bool _use_gsl_routine_for_high_x{false}; //!< flag for whether to use the gsl routine for high x
		gsl::workspace_type _workspace{nullptr}; //!< gsl workspace for calling gsl integration routines
		std::vector<gsl::workspace_type> _workspaces; //!< gsl workspaces for calling gsl integration routines

		static constexpr uint DEFAULT_GAULEG_POINTS = 50;
	public:
		Grid() = delete; //!< default constructor deleted; must provide information to fill the grid
		/**
		 *  Fills the grid with @a nx grid points according to @a grid_fill_type and sets up @a gauss_points gauss-legendre points
		 *  @param xtab Array of tabulated grid points to ensure the grid contains for easy retrieval
		 *  @param nx number of grid points
		 *  @param grid_fill_type Specifies how the grid is laid out
		 */
		Grid(grid_type const& xtab, uint nx, uint grid_fill_type=LOG_LIN);
		~Grid() = default; //!< default destructor

		explicit Grid(Grid& other); //!< constructs grid from another
		void operator=(Grid& other); //!< assigns grid from another

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

		inline bool splitIntervals() const { return _split_interval; }
		
		/** Getter for the gauss-legendre abscissae */
		inline std::vector<gauleg_type> const& abscissae() const { return _Xi; }
		/** Getter for the gauss-legendre weights */
		inline std::vector<gauleg_type> const& weights() const { return _Wi; }

		/** Getter for the grid size */
		inline uint size() const { return _points.size(); }

		/** Getter for the ntab array */
		inline ntab_type const& ntab() const { return _ntab; }
		/** Getter for the ntab array */
		inline ntab_type& ntab() { return _ntab; }

	
		/** Const iterator to the beginning of the underlying array */
		inline grid_type::const_iterator begin() const { return _points.begin(); }
		/** Iterator to the beginning of the underlying array */
		inline grid_type::iterator begin() { return _points.begin(); }
		/** Const iterator to the end of the underlying array */
		inline grid_type::const_iterator end() const { return _points.end(); }
		/** Iterator to the end of the underlying array */
		inline grid_type::iterator end() { return _points.end(); }

		/**
		 *  @brief Splits the grid intervals into the provided intervals.
		 *  @param intervals A list of points at which to split the interval
		 *  @param sizes Sizes of each interval
		 */
		void splitConvolution(
			std::vector<double> const& intervals,
			std::vector<double> const& sizes);

		/** @brief sets a flag that uses a higher accuracy (but slower) GSL routine for x>0.8 */
		inline void useGSLRoutineForHighX() { _use_gsl_routine_for_high_x = true; }

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
		void initGridLog(grid_type const& xtab, uint nx);
		/** Fills the grid with log-spacing and a calculated set of additional linear points from \f$0.1<x<1.0\f$ */
		void initGridLogLin(grid_type const& xtab, uint nx);
		
		/** Initializes the set of gauss-legendre weights and abscissae */
		void initGauLeg(double x1, double x2, std::vector<double> & Xi, std::vector<double> & Wi);
	};
}

#endif // __GRID_HPP
