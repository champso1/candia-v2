/**
 *  @file ArrayGrid.hpp
 *  @brief Contains the @a ArrayGrid class which is nothing more than an array with some convenience methods
 */

#pragma once

#include <initializer_list>

#include "Candia-v2/Common.hpp"

namespace Candia2
{
	/**
	 *  @brief stores an array and provides some convenience methods
	 */
	class ArrayGrid final
	{
	public:
		using value_type = double;
		using size_type = std::size_t; //!< alias for a size type
		using base_type = std::vector<value_type>; //!< alias for the underlying array type

	private:
		base_type _base{}; //!< underlying array
	public:
		/** @brief default constructor */
		ArrayGrid() {}
		/**
		 *  @brief default Initializes the array with zeros
		 *  @param size Size of the array
		 */
		explicit ArrayGrid(size_type size) : _base(size, 0.0) {}
		/**
		 *  @brief default Initializes the array with a set of points
		 *  @param points Another set of points
		 */
		explicit ArrayGrid(base_type const& points) : _base{points} {}
		/**
		 *  @brief default initializes the array with another @a ArrayGrid
		 *  @param other The other @a ArrayGrid
		 */
		explicit ArrayGrid(ArrayGrid const& other) : _base(other._base) {}
		explicit ArrayGrid(ArrayGrid&& other) = default;

		/**
		 *  @brief Initializes using list-initialization
		 *  @param l an initializer list of values to initialize the array with
		 */
		ArrayGrid(std::initializer_list<value_type> l) : _base(l) {}

		/** Copy assignment operator performs like the copy constructor */
		inline void operator=(ArrayGrid const& other) { _base = other.base(); }

		/** @brief returns the size of the array */
		inline size_type size() const noexcept { return _base.size(); }

		/** getter for the underlying array */
		inline base_type  const& base() const noexcept { return _base; }
		/** Zeros the entire array and performs some cleanup. */
		inline void zero() noexcept { for (double& x : _base) x = 0; }
		
		inline value_type operator[](size_type idx) const { return _base[idx]; }
		inline value_type& operator[](size_type idx) { return _base[idx]; }

		/** const begin iterator */
		inline base_type::const_iterator begin() const { return _base.cbegin(); }
		/** const end iterator */
		inline base_type::const_iterator end() const { return _base.cend(); }
		/** begin iterator */
		inline base_type::iterator begin() { return _base.begin(); }
		/** end iterator */
		inline base_type::iterator end() { return _base.end(); }
	};

	/**
	 *  @brief Templated class to ease in specifying nested @a ArrayGrid
	 */
	template <uint N>
	struct MultiDimArrayGrid
	{
		typedef typename MultiDimArrayGrid<N-1>::type Nested;
		typedef std::vector<Nested> type;
	};
	template<>
	struct MultiDimArrayGrid<1>
	{
		typedef std::vector<ArrayGrid> type;
	};

    template <uint N>
	using MultiDimArrayGrid_t = MultiDimArrayGrid<N>::type;


	/**
	 *  @brief second iteration on the "arraygrid"
	 *  this second version uses only 1 backing array for all dimensions,
	 *  rather than separate allocations for all
	 *  and the accessor method (operator()) handles strides of various dimensions
	 *  this should severely reduce cache-misses/page-faults for accessing adjacent elements,
	 *  which should make things much faster
	 *  the only negative would be that now a larger amount of contiguous memory is required,
	 *  but this is only a few hundred MBs at N3LO for the highest precision,
	 *  and if you're planning on doing that you would be doing it on a suitable machine anyway
	 *  worst case the operating system would just complain, it shouldn't explode your machine
	 */
	template <typename T, uint D>
	class ArrayGrid2Base final
	{
	public:
		using value_type = T;
	private:
		std::vector<value_type> _data{};
		std::array<uint, D> _sizes{}, _strides{};
	public:
		ArrayGrid2Base(std::array<uint, D> const& sizes)
			: _sizes{sizes}
		{
			uint current_stride = 1;
			for (int i=D-1; i>=0; --i) {
				_strides[i] = current_stride;
				current_stride *= sizes[i];
			}
			_data.resize(current_stride);
			std::ranges::fill(_data, double{0});
		}

		template <typename... TArgs, uint... ints>
		uint acquire_idx([[maybe_unused]] std::integer_sequence<uint, ints...> intseq, TArgs... args) {
			return ((_strides[ints]*args) + ...);
		}

		inline auto size(uint idx) const { return _sizes[idx]; }
		inline auto const& sizes() const { return _sizes; }

		inline void resize(std::array<uint, D> const& new_sizes) {
			_data.clear();
			uint current_stride = 1;
			for (int i=D-1; i>=0; --i) {
				_strides[i] = current_stride;
				current_stride *= new_sizes[i];
			}
			_data.resize(current_stride);
			std::ranges::fill(_data, double{0});
			std::ranges::copy(new_sizes, _sizes.begin());
		}
		
		template <typename... TArgs>
		decltype(auto) operator()(TArgs... args)
		{
			constexpr auto num_args = sizeof...(args);
			uint idx = acquire_idx(std::make_integer_sequence<uint, num_args>{}, std::forward<TArgs>(args)...);
			if constexpr (sizeof...(args) == D-1)
				return std::span<double>(_data.begin() + idx, _sizes.back());
			else if constexpr (sizeof...(args) == D)
				return _data[idx];
			throw std::runtime_error("cannot access anything other than last dimension (span) or singular element");
		}

		template <typename... TArgs>
		void copy_from(std::span<double> other, TArgs... args)
		{
			static_assert(sizeof...(args) == D-1, "must copy arraygrid to final dimension (i.e. provide D-1 arguments)");
			if (other.size() != _sizes.back())
				log(LOG_ERROR, "arraygrid2", "size mismatch during copy ({} vs. {})", other.size(), _sizes.back());
		    auto this_span = this->operator()(std::forward<TArgs>(args)...);
			std::ranges::copy(other, this_span.begin());
		}
	};

	template <uint D>
	using ArrayGrid2 = ArrayGrid2Base<double, D>;
}
