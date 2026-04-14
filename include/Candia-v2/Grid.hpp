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

		static inline int error_print_count = 35;

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

	/** @brief Struct that provides the options for the behavior of the @a Grid */
	struct GridOptions final
	{
		bool use_gsl_conv_routine{false};  //!< flag for whether to use a gsl routine to perform the convolution
		bool use_gsl_interp_routine{true}; //!< flag for whether to use a gsl routine to perform the interpolation
		bool use_alt_mapping{false}; //!< flag for whether to use a different mapping for convolutions to try and increase accuracy
	};

	/** @brief Base class for types of grid fillers. */
	struct GridFillerBase
	{
		std::vector<double> xtab{};
		std::vector<int> ntab{};
		double min{1e-5};

		GridFillerBase() = default;
		GridFillerBase(double min_)
			: min{min_}
		{}
		virtual ~GridFillerBase() = default;
		
		virtual uint fill(std::vector<double>& points) = 0;
	};

	/** @brief Fills a grid linearly with a particular number of points. */
	struct GridFillerLin final : public GridFillerBase
	{
		uint size{601};
		
		GridFillerLin() = default;
		GridFillerLin(double min, uint size_)
			: GridFillerBase(min),
			  size{size_} {}
		uint fill(std::vector<double>& points) override;
	};

	/** @brief Fills a grid logarithmically with a particular number of points. */
	struct GridFillerLog final : public GridFillerBase
	{
		uint size{501};
		
		GridFillerLog() = default;
		GridFillerLog(double min, uint size_)
			: GridFillerBase(min),
			  size{size_} {}
		uint fill(std::vector<double>& points) override;
	};

	/** @brief Fills a grid logarithmically in the lower interval and linearly in the upper interval. */
	struct GridFillerLogLin final : public GridFillerBase
	{
		uint log_size{151};
		uint lin_size{151};
		double pivot{0.1};

		GridFillerLogLin() = default;
		GridFillerLogLin(double min, uint log_size_, uint lin_size_, double pivot_)
			: GridFillerBase(min),
			  log_size{log_size_}, lin_size{lin_size_}, pivot{pivot_} {}
		uint fill(std::vector<double>& points) override;
	};

	/**
	 *  @brief Fills a grid logarithmically in the lower interval, linearly in a middle interval,
	 *  and quadratically in the upper interval, packing points at x->1.
	 */
	
	struct GridFillerLogLinQuad final : public GridFillerBase
	{
		uint log_size{101};
		uint lin_size{51};
		uint quad_size{26};
		double pivot1{0.1}, pivot2{0.9};

		GridFillerLogLinQuad() = default;
		GridFillerLogLinQuad(double min, uint log_size_, uint lin_size_, uint quad_size_)
			: GridFillerBase(min),
			  log_size{log_size_}, lin_size{lin_size_}, quad_size{quad_size_} {}
		uint fill(std::vector<double>& points) override;
	};

	/** @brief helper function to create a grid filler object for the @a Grid constructor */
    template <typename TGridFiller, typename... TGridFillerArgs>
	inline std::unique_ptr<GridFillerBase> make_grid_filler(TGridFillerArgs&&... args)
	{
		return std::make_unique<TGridFiller>(std::forward<TGridFillerArgs>(args)...);
	}


	struct GaussLegendreArgs final
	{
		uint default_gauss_points{100};
		bool split_interval{false};
		std::vector<std::tuple<double,double,uint>> intervals{};
	};

	/**
	 *  @brief Class that contains the interpolation/convolution grid and methods
	 */
	class Grid final : public OptionsBase<GridOptions>
	{
    public:
		using value_type = double; //!< alias for underlying grid type. not settable at this point
		using grid_type = std::vector<value_type>; //!< alias for the underlying grid type
		using gauleg_type = std::vector<value_type>; //!< alias for the type of the array of gauss-legendre weights/abscissae
		using ntab_type = std::vector<int>; //!< alias for the type of the calulated ntab array
		using grid_filler_type = std::unique_ptr<GridFillerBase>; //!< alias for the type that fills the grid
		using gsl_conv_errors = std::vector<std::tuple<value_type,value_type,value_type>>; //!< alias for list of errors obtained after gsl integration
		using mapping_type = std::vector<std::function<std::pair<double,double>(double,double)>>; //!< type for mappings

		/** @brief struct passed to gsl integration routine */
		struct GSLIntegrationParams final
		{
			Grid& g;

			value_type x;
			uint k;
			value_type logx;
			value_type eplus1;

			ArrayGrid& A;
			Expression& E;
		};

	private:
		grid_type _points{}; //!< grid points
		ntab_type _ntab{};     //!< stored indices for the tabulated grid points
		grid_type _xtab{};     //!< stored values of the tabulated grid points

		GaussLegendreArgs _gauleg_args{}; //!< arguments/inputs/setup for gauleg integration
		std::vector<gauleg_type> _Xi{}; //!< list of split-up gauleg abscissae per interval
		std::vector<gauleg_type> _Wi{}; //!< list of split-up gauleg weights per interval

		std::vector<gsl::workspace_type> _workspaces; //!< gsl workspaces for calling gsl integration routines
		std::vector<gsl::interp_type> _interps; //!< gsl interp objects for interpolation
		std::vector<gsl::interp_accel_type> _interp_accels; //!< interp interpolation acceleration objects
		gsl_conv_errors _gsl_conv_errors{}; //!< stored values of GSL fails (if any); for debugging purposes

		std::vector<std::function<std::pair<double,double>(double,double)>> _mappings{}; //!< mapping functions

		bool _use_cached_expressions{false}; //!< whether to use cached expressions. set via DGLAPSolver based on its options

	public:
		struct EnumerateIterator final
		{
			grid_type const& data;

			struct Iterator final
			{
				grid_type const& v;
				grid_type::size_type idx;

				inline auto operator*() const {
					return std::pair<grid_type::size_type, grid_type::value_type>{idx, v[idx]};
				}

				inline Iterator& operator++() { ++idx; return *this; }
				inline bool operator!=(Iterator const& other) { return idx != other.idx; }
			};

			inline auto begin() { return Iterator{data, 0}; }
			inline auto end() { return Iterator{data, data.size()}; }
		};

		inline auto enumerate() const {
			return EnumerateIterator{_points};
		}
		
	public:
		Grid() = delete; //!< default constructor deleted; must provide information to fill the grid
		/**
		 *  fills the grid, sets up GSL objects, and other initialization functions
		 *  @param xtab Array of tabulated grid points to ensure the grid contains for easy retrieval
		 *  @param grid_filler an object that will fill the grid in a particular way
		 *  @param gauleg_args arguments to setup/initialize how gauleg integration behaves
		 */
		Grid(grid_type const& xtab, grid_filler_type grid_filler, GaussLegendreArgs const& gauleg_args);
		~Grid() = default;

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
			Expression& E, ArrayGrid& A,
			value_type eplus1,
			gauleg_type const& X, gauleg_type const& W);

		/** @brief sets the flag to use the cached expressions */
		void useCachedExpressions() { _use_cached_expressions = true; }

		/** @brief returns the flag for whether to use the cached expressions */
		bool usingCachedExpressions() { return _use_cached_expressions; }
		
		/** @brief returns the registered mapping types */
		mapping_type const& getMappings();

		/** @brief called from DGLAPSolver::evolve() to ensure that the mappings are setup correctly */
		void setupMappings() { UNUSED(getMappings()); }

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
		value_type interpolate(ArrayGrid& y, value_type x);
		/**
		 *  @brief Performs a convolution between an array @a A and an expression @a E
		 *  @param A the array
		 *  @param E the expression
		 *  @param k the grid index to perform the convolution at
		 */
		value_type convolution(ArrayGrid& A, Expression& E, uint k);

		/** @brief retrieves info on all GSL errors, if any */
		inline auto const& getGSLConvolutionErrors() const { return _gsl_conv_errors; }
	private:
	    /** Initializes the set of gauss-legendre weights and abscissae */
		void initGauLeg(value_type x1, value_type x2, gauleg_type& Xi, gauleg_type& Wi);
	};
}

#endif // __GRID_HPP
