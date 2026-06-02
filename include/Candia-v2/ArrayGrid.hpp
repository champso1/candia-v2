/**
 *  @file ArrayGrid.hpp
 *  @brief Contains the @a ArrayGrid class which is nothing more than an array with some convenience methods
 */

#pragma once

#include <vector>
#include <array>

#include "Candia-v2/Common.hpp"

namespace Candia2
{
	/**
	 *  @brief second iteration on the "arraygrid"
	 *  this second version uses only 1 backing array for all dimensions,
	 *  rather than separate allocations for all
	 *  and the accessor method (operator()) handles strides of various dimensions
	 *  this should severely reduce cache-misses/page-faults for accessing adjacent elements,
	 *  the only negative would be that now a larger amount of contiguous memory is required,
	 *  but this is only a few hundred MBs at N3LO for the highest precision,
	 *  and if you're planning on doing that you would be doing it on a suitable machine anyway
	 *  worst case the operating system would just complain, it shouldn't explode your machine
	 */
	template <typename T, uint D>
	class ArrayGridBase final
	{
	public:
		using value_type = T;
	    using data_type = std::vector<value_type>;
		using size_type = data_type::size_type;
	private:
		std::vector<value_type> _data{};
		std::array<uint, D> _sizes{}, _strides{};
	public:
		ArrayGridBase() = default;
		ArrayGridBase(std::array<uint, D> const& sizes)
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
		uint acquire_idx([[maybe_unused]] std::integer_sequence<uint, ints...> intseq, TArgs... args)
		{
			return ((_strides[ints]*args) + ...);
		}

		inline auto size(uint idx) const { return _sizes[idx]; }
		inline auto const& sizes() const { return _sizes; }

		inline void resize(std::array<uint, D> const& new_sizes)
		{
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

		inline void clear()
		{
			_data.clear();
			std::ranges::fill(_sizes, 0.0);
			std::ranges::fill(_strides, 0.0);
		}

		auto begin() requires (D==1) { return _data.begin(); }
		auto end()   requires (D==1) { return _data.end(); }
		auto begin() const requires (D==1) { return _data.begin(); }
		auto end()   const requires (D==1) { return _data.end(); }
		
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
	};

	template <typename T>
	class ArrayGridBase<T, 1> final
	{
	public:
		using value_type = T;
	private:
		std::vector<value_type> _data{};
	public:
		ArrayGridBase() = default;
		virtual ~ArrayGridBase() = default;
		explicit ArrayGridBase(uint size) : _data(size, 0.0) {}
		ArrayGridBase(ArrayGridBase const& other)
			: _data{other._data}
		{}
		ArrayGridBase(ArrayGridBase&& other)
			: _data{other._data}
		{
			other._data.clear();
		}
		template <typename TIterator>
		ArrayGridBase(TIterator begin, TIterator end)
		{
			_data = std::vector<value_type>(begin, end);
		}

		template <typename TContainer>
		ArrayGridBase(TContainer const& container)
			: _data(container)
		{}

		ArrayGridBase(std::span<value_type> view)
			: _data(view.begin(), view.end())
		{}

		void operator=(ArrayGridBase const& other)
		{
			_data.resize(other._data.size());
			std::ranges::copy(other._data, _data.begin());
		}
		void operator=(ArrayGridBase&& other)
		{
			_data.resize(other._data.size());
			std::ranges::copy(other._data, _data.begin());
			other._data.clear();
		}

		inline auto size() const { return _data.size(); }
		inline auto empty() const { return _data.empty(); }

		inline void resize(uint new_size)
		{
			_data.clear();
			_data.resize(new_size);
			std::ranges::fill(_data, double{0});
		}

		inline void clear()
		{
			_data.clear();
		}

		auto begin() { return _data.begin(); }
		auto end()   { return _data.end(); }
		auto begin() const { return _data.begin(); }
		auto end()   const { return _data.end(); }
		
	    value_type operator()(uint idx) const { return _data[idx]; }
		value_type operator[](uint idx) const { return _data[idx]; }
		value_type& operator()(uint idx) { return _data[idx]; }
		value_type& operator[](uint idx) { return _data[idx]; }

		auto const& data() const { return _data; }
		auto& data() { return _data; }
		std::span<value_type> view() { return std::span<double>(_data); }
	};
	

	template <uint D>
	using ArrayGridN = ArrayGridBase<double, D>;
	using ArrayGrid = ArrayGridBase<double, 1>;
    using ArrayGridView = std::span<double>;
}
