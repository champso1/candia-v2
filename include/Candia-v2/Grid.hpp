/**
 *  @file Grid.hpp
 *  @brief Contains the @a Grid class which contains the interpolation/convolution grids and associated routines.
 */

#pragma once

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"
#include "Candia-v2/ArrayGrid.hpp"

#include <cmath>

namespace Candia2
{
    /** @brief Base class for types of grid fillers. */
	struct GridFillerBase
	{
		using mapping_function_type = std::function<std::pair<double,double>(double,double)>;
		
		std::vector<double> xtab{};
		std::vector<int> ntab{};
		double min{1e-5};

		GridFillerBase() = default;
		GridFillerBase(double min_)
			: min{min_}
		{}
		virtual ~GridFillerBase() = default;
		
		virtual uint fill(std::vector<double>& points) = 0;
		virtual std::span<mapping_function_type> getMappings(double x) = 0;
	};

	/** @brief Fills a grid linearly with a particular number of points. */
	struct GridFillerLin final : public GridFillerBase
	{
	private:
		static std::vector<mapping_function_type> _mappings;
		static std::span<mapping_function_type> _mapping_span;
	public:
		uint size{601};
		
		GridFillerLin() = default;
		GridFillerLin(double min, uint size_)
			: GridFillerBase(min),
			  size{size_} {}
		uint fill(std::vector<double>& points) override;

		inline std::span<mapping_function_type> getMappings([[maybe_unused]] double x) override { return _mapping_span; }
	};

	/** @brief Fills a grid logarithmically with a particular number of points. */
	struct GridFillerLog final : public GridFillerBase
	{
	private:
		static std::vector<mapping_function_type> _mappings;
		static std::span<mapping_function_type> _mapping_span;
	public:
		uint size{501};
		
		GridFillerLog() = default;
		GridFillerLog(double min, uint size_)
			: GridFillerBase(min),
			  size{size_} {}
		uint fill(std::vector<double>& points) override;
		inline std::span<mapping_function_type> getMappings([[maybe_unused]] double x) override { return _mapping_span; }
	};

	/** @brief Fills a grid logarithmically in the lower interval and linearly in the upper interval. */
	struct GridFillerLogLin final : public GridFillerBase
	{
	public:
		uint log_size{151};
		uint lin_size{151};
		double pivot{0.1};
	private:
		std::vector<mapping_function_type> _mappings;
		void setupMappings()
		{
			_mappings = {
				[&](double x, double z) {
					auto a = pivot/x;
					return std::make_pair(x*std::pow(a, z), x*std::pow(a, z)*std::log(a));
				},
				[&]([[maybe_unused]] double x, double z) { return std::make_pair(pivot+(1.0-pivot)*z, 1.0-pivot); },
				[](double x, double z) { return std::make_pair(x+(1.0-x)*z, 1.0-x); },
			};
		}
	public:
		GridFillerLogLin()
		{
			setupMappings();
		}
		GridFillerLogLin(double min, uint log_size_, uint lin_size_, double pivot_)
			: GridFillerBase(min),
			  log_size{log_size_}, lin_size{lin_size_}, pivot{pivot_}
		{
			setupMappings();
		}
		uint fill(std::vector<double>& points) override;
		
		inline std::span<mapping_function_type> getMappings(double x) override
		{
			if (x < 0)
				return std::span(_mappings.begin(), _mappings.end());
			else if (x > 0 && x < pivot)
				return std::span(_mappings.begin(), _mappings.begin()+2);
			else
				return std::span(_mappings.begin()+2, _mappings.end());
		}
	};

	/**
	 *  @brief Fills a grid logarithmically in the lower interval, linearly in a middle interval,
	 *  and quadratically in the upper interval, packing points at x->1.
	 */
	
	struct GridFillerLogLinQuad final : public GridFillerBase
	{
	public:
		uint log_size{101};
		uint lin_size{51};
		uint quad_size{26};
		double pivot1{0.1}, pivot2{0.9};
	private:
		std::vector<mapping_function_type> _mappings;
		void setupMappings()
		{
			_mappings = {
				[&](double x, double z) {
					auto a = pivot1/x;
					return std::make_pair(x*std::pow(a, z), x*std::pow(a, z)*std::log(a));},
				[&]([[maybe_unused]] double x, double z) { return std::make_pair(pivot1+(pivot2-pivot1)*z, pivot2-pivot1); },
				[&]([[maybe_unused]] double x, double z) { return std::make_pair(1.0-(1.0-pivot2)*std::pow(1.0-z, 3), 3.0*(1.0-pivot2)*(1.0-z)*(1.0-z)); },
				[&](double x, double z) { return std::make_pair(x+(pivot2-x)*z, pivot2-x); },
				[&]([[maybe_unused]] double x, double z) { return std::make_pair(1.0-(1.0-pivot2)*std::pow(1.0-z, 3), 3.0*(1.0-pivot2)*(1.0-z)*(1.0-z)); },
				[](double x, double z) { return std::make_pair(1.0-(1.0-x)*std::pow(1.0-z, 3), 3.0*(1.0-x)*(1.0-z)*(1.0-z)); },
			};
		}
	public:
		GridFillerLogLinQuad()
		{
			setupMappings();
		}
		GridFillerLogLinQuad(double min, uint log_size_, uint lin_size_, uint quad_size_)
			: GridFillerBase(min),
			  log_size{log_size_}, lin_size{lin_size_}, quad_size{quad_size_}
		{
			setupMappings();
		}
		uint fill(std::vector<double>& points) override;

		inline std::span<mapping_function_type> getMappings(double x) override
		{
			if (x < 0)
				return std::span(_mappings.begin(), _mappings.end());
			else if (x > 0 && x < pivot1)
				return std::span(_mappings.begin(), _mappings.begin()+3);
			else if (x >= pivot1 && x < pivot2)
				return std::span(_mappings.begin()+3, _mappings.begin()+3+2);
			else
				return std::span(_mappings.begin()+3+2, _mappings.end());
		}
	};

	/** @brief helper function to create a grid filler object for the @a Grid constructor */
    template <typename TGridFiller, typename... TGridFillerArgs>
	inline std::unique_ptr<GridFillerBase> make_grid_filler(TGridFillerArgs&&... args)
	{
		return std::make_unique<TGridFiller>(std::forward<TGridFillerArgs>(args)...);
	}



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
		using value_type = double; //!< alias for underlying grid type. not settable at this point
		using grid_type = std::vector<value_type>; //!< alias for the underlying grid type
		using gauleg_type = std::vector<value_type>; //!< alias for the type of the array of gauss-legendre weights/abscissae
		using ntab_type = std::vector<int>; //!< alias for the type of the calulated ntab array
	private:
		std::reference_wrapper<GridFillerBase> _filler;
		grid_type _points{}; //!< grid points
		ntab_type _ntab{};     //!< stored indices for the tabulated grid points
		grid_type _xtab{};     //!< stored values of the tabulated grid points

		ConvIntArgs _convint_args; //!< contains misc convolution/interpolation options/args
		gauleg_type _Xi{}; //!< list of split-up gauleg abscissae per interval
		gauleg_type _Wi{}; //!< list of split-up gauleg weights per interval

		bool _use_cached_expressions{false}; //!< whether to use cached expressions. set via DGLAPSolver based on its options
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
		Grid() = delete; //!< must provide information to fill the grid
		/**
		 *  fills the grid, sets up GSL objects, and other initialization functions
		 *  @param xtab Array of tabulated grid points to ensure the grid contains for easy retrieval
		 *  @param grid_filler an object that will fill the grid in a particular way
		 *  @param gauleg_args arguments to setup/initialize how gauleg integration behaves
		 */
		Grid(grid_type const& xtab, GridFillerBase& grid_filler, ConvIntArgs const& gauleg_args);
		Grid(Grid const& other) = default;
		Grid(Grid&& other) = default;
		~Grid() = default;

		inline uint size() const { return _points.size(); }

		inline GridFillerBase& filler() { return _filler.get(); }

		inline grid_type& xtab() { return _xtab; }
		inline grid_type const& xtab() const { return _xtab; }
		inline ntab_type const& ntab() const { return _ntab; }
		inline ntab_type& ntab() { return _ntab; }

		inline grid_type const& points() const { return _points; }
		inline value_type at(uint idx) const { return _points.at(idx); };
		inline value_type operator[](uint idx) const { return _points[idx]; }

		inline gauleg_type const& abscissae() const { return _Xi; }
		inline gauleg_type const& weights() const { return _Wi; }

		inline grid_type::const_iterator begin() const { return _points.begin(); }
		inline grid_type::iterator begin() { return _points.begin(); }
		inline grid_type::const_iterator end() const { return _points.end(); }
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
	private:
	    /** Initializes the set of gauss-legendre weights and abscissae */
		void initGauLeg(value_type x1, value_type x2, gauleg_type& Xi, gauleg_type& Wi);
	};
}
