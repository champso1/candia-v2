#ifndef __FUNCARRGRID_HPP
#define __FUNCARRGRID_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Grid.hpp"

#include <functional>
#include <memory>
#include <utility>
#include <map>

namespace Candia2
{
	class ArrayGrid final
	{
	private:
		std::reference_wrapper<const Grid> _grid;

		using base_type = std::vector<double>;
		using cache_type = std::map<double, double>;
		using size_type = std::size_t;

		base_type _base{};
		cache_type _cache{};
	public:
		ArrayGrid(Grid const& grid)
			: _grid{std::cref(grid)}, _base(grid.size(), 0.0) {}
		ArrayGrid(Grid const& grid, base_type const& points)
			: _grid{std::cref(grid)}, _base{points} {}
		~ArrayGrid() = default;

		inline base_type  const& base() const noexcept { return _base; }
		inline cache_type const& cache() const noexcept { return _cache; }
	    size_type size() const noexcept;
		void zero() noexcept;
		
		double& operator[](uint idx); // base accessor for points on the grid
		double& operator()(double x); // accessor for points off the grid

	private:
		void addPoint(double x);
		double interpolate(double x);
	};

	

	class FunctionGrid final
	{
	private:
		using cache_type = std::map<double, std::array<double, 3>>;
		using grid_type = std::reference_wrapper<const Grid>;
		using expr_type = std::unique_ptr<Expression>;

		cache_type _cache{}; //!< stores evaluated pieces of the function
		grid_type _grid; //!< underlying grid
		expr_type _func; //!< unique_ptr to corresponding function
	public:
		FunctionGrid(Grid const& grid, expr_type expr);
		FunctionGrid(FunctionGrid&& other)
			: _cache{std::move(other._cache)},
			  _grid{std::move(other._grid)},
			  _func{std::move(other._func)}
		{ }
		~FunctionGrid() = default;

		FunctionGrid(FunctionGrid const& other) = delete;
		void operator=(FunctionGrid const& other) = delete;

		enum FunctionPart : uint
		{
			REGULAR = 0,
			PLUS = 1,
			DELTA = 2
		};

		double operator()(double x, uint function_part);
		void addFunctionPoints(std::vector<double> const& X);
		void addFunctionPoint(double x);

		double convolution(ArrayGrid & A, uint k, bool split_interval=false);

		inline cache_type const& cache() const noexcept { return _cache; }
	};



	template <uint N>
	struct MultiDimArrayGrid
	{
		typedef typename MultiDimArrayGrid<N-1>::type Nested;
		typedef std::vector<Nested> type;
	};
	template <>
	struct MultiDimArrayGrid<1>
	{
		typedef std::vector<ArrayGrid> type;
	};

	template <uint N>
	using MultiDimArrayGrid_t = MultiDimArrayGrid<N>::type;
}

#endif // __FUNCARRGRID_HPP
