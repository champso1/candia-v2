#ifndef __FUNCARRGRID_HPP
#define __FUNCARRGRID_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Grid.hpp"

#include <functional>
#include <fstream>
#include <utility>
#include <map>

namespace Candia2
{
	class ArrayGrid final
	{
	private:
		// underlying grid object
		std::reference_wrapper<const Grid> _grid;

		// we keep the set of existing points,
		// which is just an array that matches with the grid
		// we need this becuase we need the interpolation routine
		// to use only the existing points
		// we shouldn't interpolate using points that were themselves interpolated
		// we store those extra points in a map
		using base_type = std::vector<double>;
		using cache_type = std::map<double, double>;
		using size_type = std::size_t;

		using const_base_it = base_type::const_iterator;
		using base_it = base_type::iterator;
		using const_cache_it = cache_type::const_iterator;
		using cache_it = cache_type::iterator;

		base_type _base{};
		cache_type _cache{};
	public:
		ArrayGrid(Grid const& grid)
			: _grid{std::cref(grid)}, _base(grid.Size(), 0.0) {}
		ArrayGrid(Grid const& grid, base_type const& points)
			: _grid{std::cref(grid)}, _base{points} {}
		~ArrayGrid() = default;

		inline base_type  const& base() const { return _base; }
		inline cache_type const& cache() const { return _cache; }
	    size_type size() const noexcept;
		void zero() noexcept;
		
		double& operator[](uint idx); // base accessor for points on the grid
		double& operator()(double x); // accessor for points off the grid

		// ==============================
		// c++ range begin/end range functions
		// ==============================
		inline const_base_it  base_begin() const { return _base.begin(); }
		inline base_it        base_begin() { return _base.begin(); }
		inline const_base_it  base_end()   const { return _base.end(); }
		inline base_it        base_end()   { return _base.end(); }
		inline const_cache_it cache_begin() const { return _cache.begin(); }
		inline cache_it       cache_begin() { return _cache.begin(); }
		inline const_cache_it cache_end()   const { return _cache.end(); }
		inline cache_it       cache_end()   { return _cache.end(); }

	private:
		void addPoint(double x);
		double interpolate(double x);
	};

	

	class FunctionGrid final
	{
	private:
		// aliases
		using cache_type = std::map<double, std::array<double, 3>>;
		using iterator = cache_type::iterator;
		using const_iterator = cache_type::const_iterator;

		// map that stores computed values of the three parts of a function,
		// that being the regular, plus, and delta pieces
		// the key is the index into the underlying grid,
		// and the value is an array of three values
		// corresponding to the
		cache_type _cache{};

		// underlying grid object
		std::reference_wrapper<const Grid> _grid;

		// alias to the used type
		std::unique_ptr<Expression> _func;
	public:
		// fills the grid with values of the function at every grid point
		// defined in the underlying grid
		FunctionGrid(Grid const& grid, std::unique_ptr<Expression> expr);
		FunctionGrid(FunctionGrid const& other) = delete;
		void operator=(FunctionGrid const& other) = delete;
		FunctionGrid(FunctionGrid&& other)
			: _cache{std::move(other._cache)},
			  _grid{std::move(other._grid)},
			  _func{std::move(other._func)}
		{ }
		~FunctionGrid() = default;

		// these are directly the indicies used to interface with the array in the map
		enum FunctionPart : uint
		{
			REGULAR = 0,
			PLUS = 1,
			DELTA = 2
		};

		inline iterator begin() { return _cache.begin(); }
		inline const_iterator begin() const { return _cache.begin(); }
		inline iterator end() { return _cache.end(); }
		inline const_iterator end() const { return _cache.end(); }

		double operator()(double x, uint function_part);
		void addFunctionPoints(std::vector<double> const& X);
		void addFunctionPoint(double x);

		double convolution(ArrayGrid & A, uint k);

		inline
		cache_type const& cache() const { return _cache; }
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
