/**
 *  @file Grid.hpp
 *  @brief Contains the @a Grid class which contains the interpolation/convolution grids and associated routines.
 */

#ifndef __GRID_HPP
#define __GRID_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"
#include "Candia-v2/Options.hpp"
#include "Candia-v2/ArrayGrid.hpp"

#include <vector>
#include <memory>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

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

		inline auto interp_deleter = [](gsl_interp* interp){ gsl_interp_free(interp); };
		using interp_deleter_type = decltype(interp_deleter);
		using interp_type = std::unique_ptr<gsl_interp, interp_deleter_type>;
		inline auto interp_accel_deleter = [](gsl_interp_accel* accel){ gsl_interp_accel_alloc(); };
		using interp_accel_deleter_type = decltype(interp_accel_deleter);
		using interp_accel_type = std::unique_ptr<gsl_interp_accel, interp_accel_deleter_type>;
		inline auto const* interp_cspline_type = gsl_interp_cspline;

		inline auto make_interp = [](uint size){ return interp_type(gsl_interp_alloc(interp_cspline_type, size), interp_deleter); };
		inline auto make_interp_accel = [](){ return interp_accel_type(gsl_interp_accel_alloc(), interp_accel_deleter); };

		static inline int error_print_count = 50;

		extern "C" {
			static inline void void_error_handler(
				const char * reason, const char * file,
				int line, int gsl_errno)
			{
				UNUSED(reason);
				UNUSED(file);
				UNUSED(line);
				UNUSED(gsl_errno);
			    return;
			}
		}

		static bool error_handler_set = (([](){ /* gsl_set_error_handler(void_error_handler); */ return; })(), true);
	}
	
	struct GridOptions final
	{
		bool use_gsl_conv_routine{false};  //!< flag for whether to use a gsl routine to perform the convolution
		bool use_gsl_interp_routine{true}; //!< flag for whether to use a gsl routine to perform the interpolation
		bool use_alt_mapping{false}; //!< flag for whether to use a different mapping for convolutions to try and increase accuracy
		
	};

	/**
	 *  @brief Class that contains the interpolation/convolution grid and the methods to perform the interpolation and convolution.
	 */
	class Grid final : public OptionsBase<GridOptions>
	{
	public:
		using value_type = double;
		using grid_type = std::vector<value_type>; //!< alias for the underlying grid type
		using gauleg_type = std::vector<value_type>; //!< alias for the type of the array of gauss-legendre weights/abscissae
		using ntab_type = std::vector<int>; //!< alias for the type of the calulated ntab array
		using arraygrid_type = ArrayGrid;
		using expression_type = Expression;

		enum GridFillType : uint
		{
			LOG = 0, //!< simple logarithmic intervals
			LOG_LIN, //!< logarithmic intervals until 0.1, then linear until 1.0
			LIN,     //!< linear all throughout
			LOG_LIN_QUAD //!< log from min-0.1, lin from 0.1-0.9, quadratic (packed at higher x) from 0.9-1.0
		};

		struct GSLIntegrationParams final
		{
			Grid& g;

			value_type x;
			uint k;
			value_type logx;
			value_type eplus1;

			arraygrid_type& A;
			Expression& E;
		};

	private:
		grid_type _points{}; //!< grid points
		ntab_type _ntab;     //!< stored indices for the tabulated grid points
		grid_type _xtab;     //!< stored values of the tabulated grid points
		using gsl_conv_errors = std::vector<std::tuple<value_type,value_type,value_type>>;

		bool _split_interval{false}; //!< flag that declares if the user has split the convolution into intervals
		std::vector<gauleg_type> _Xi{}; //!< list of split-up gauleg abscissae per interval
		std::vector<gauleg_type> _Wi{}; //!< list of split-up gauleg weights per interval
		std::vector<uint> _interval_sizes{}; //!< number of points per interval
		std::vector<gsl::workspace_type> _workspaces; //!< gsl workspaces for calling gsl integration routines
		std::vector<gsl::interp_type> _interps; //!< gsl interp objects for interpolation
		std::vector<gsl::interp_accel_type> _interp_accels; //!< interp interpolation acceleration objects 

		static constexpr uint DEFAULT_GAULEG_POINTS = 100; //!< default number of gauss-legendre points to place in the interval
		gsl_conv_errors _gsl_conv_errors{}; //!< stored values of information whenever GSL fails to perform the integration to the requested accuracy
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
		inline value_type at(uint idx) const { return _points.at(idx); };
		/** Getter for point on the grid */
		inline value_type operator[](uint idx) const { return _points[idx]; }

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
		 *
		 *  These may be left empty to not split the intervals, in which case
		 *  one single interval will be used.
		 */
		void splitConvolution(
			std::vector<value_type> const& intervals={},
			std::vector<value_type> const& sizes={});

		/**
		 *  @brief Handles a convolution with simple mappings for y -> z
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
			uint k, value_type x, auto&& yandjaccessor,
			expression_type& E, arraygrid_type& A,
			value_type eplus1,
			gauleg_type const& X, gauleg_type const& W);

		/**
		 *  @brief Uses a binary search to find the grid point closest to the given value of x
		 *  @param x value to search for
		 */
		int interpFindIdx(value_type x);
		/**
		 *  @brief Interpolates array @a y at @a x on the grid.
		 *  @param y the array grid to interpolate
		 *  @param x the value of x to interpolate at
		 */
		value_type interpolate(arraygrid_type& y, value_type x);
		/**
		 *  @brief Performs a convolution between an array @a A and an expression @a E
		 *  @param A the array
		 *  @param E the expression
		 *  @param k the grid index to perform the convolution at
		 */
		value_type convolution(arraygrid_type& A, expression_type& E, uint k);

		/** @brief retrieves info on all GSL errors, if any */
		inline auto const& getGSLConvolutionErrors() const { return _gsl_conv_errors; }
	private:
		/** Fills the grid according to the original candia-v2 method (log-spaced) */
		void initGridLog(grid_type const& xtab, uint nx);
		/** Fills the grid with log-spacing and a calculated set of additional linear points from \f$0.1<x<1.0\f$ */
		void initGridLogLin(grid_type const& xtab, uint nx);
		/** Fills the grid with linear spacing */
		void initGridLin(grid_type const& xtab, uint nx);
		/** Fills the grid with part logarithmic, part linear, and part quadratic */
		void initGridLogLinQuad(grid_type const& xtab, uint nx);
		
		/** Initializes the set of gauss-legendre weights and abscissae */
		void initGauLeg(value_type x1, value_type x2, gauleg_type& Xi, gauleg_type& Wi);
	};
}

#endif // __GRID_HPP
