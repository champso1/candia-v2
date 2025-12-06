#ifndef __COMMON_HPP
#define __COMMON_HPP

#include <vector>
#include <concepts>
#include <ranges>
#include <numbers>

typedef unsigned uint;


#define UNUSED(x) (void)(x)

#define MAX(x,y) ((x)>=(y)) ? (x) : (y)
#define MIN(x,y) ((x)<=(y)) ? (x) : (y)

#define EPS 1e-8


/// @brief 
namespace Candia2
{
	constexpr const double CF    = 4.0/3.0;
	constexpr const double NC    = 3.0;
	constexpr const double TR    = 0.5;

	constexpr const double MZ    = 91.1876;

	constexpr const double PI = std::numbers::pi;
	constexpr const double PI_2 = PI*PI;
	constexpr const double PI_3 = PI*PI*PI;

	constexpr const double Zeta2 = PI*PI/6.0;
	constexpr const double Zeta3 = 1.2020569031595942854;

	constexpr const uint DISTS = 37;
	constexpr const uint GAUSS_POINTS = 30;
	constexpr const uint INTERP_POINTS = 6;
	constexpr const uint DEFAULT_ITERATIONS = 10;
	constexpr const uint DEFAULT_TRUNC_IDX = 5;



	template <typename T, uint N>
	struct MultiDimVector
	{
		typedef typename MultiDimVector<T,N-1>::type Nested;
		typedef std::vector<Nested> type;
	};
	template <typename T>
	struct MultiDimVector<T,1>
	{
		typedef std::vector<T> type;
	};


	template<typename T>
	concept Iterable = requires(T t) {
		{ std::ranges::begin(t) } -> std::input_or_output_iterator;
		{ std::ranges::end(t) } -> std::sentinel_for<decltype(std::ranges::begin(t))>;
	};

	template <typename T>
	concept Arithmetic = std::integral<T> or std::floating_point<T>;
};


#endif // __COMMON_HPP
