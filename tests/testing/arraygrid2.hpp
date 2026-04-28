#pragma once

#include "Candia-v2/Common.hpp"
#include <sstream>
#include <utility>

namespace Candia2
{
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
