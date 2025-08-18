/** @file
 *
 *  Provides some common utilities/general functionality
 *  for the rest of the program.
 */

#ifndef __COMMON_HPP
#define __COMMON_HPP

#include <cmath>
#include <vector>
#include <concepts>
#include <ranges>
#include <numbers>

typedef unsigned int uint;


#define UNUSED(x) (void)(x)

#define MAX(x,y) ((x)>=(y)) ? (x) : (y)
#define MIN(x,y) ((x)<=(y)) ? (x) : (y)


/// @brief 
namespace Candia2
{

	/** @name Group-theoretical constants.
	 */
	///@{
	const double CF    = 4.0/3.0;
	const double NC    = 3.0;
	const double TR    = 0.5;
	///@}


	/** @brief Mass of Z-boson */
	const double MZ    = 91.1876;

	// some shorthand
	constexpr const double PI = std::numbers::pi;
	constexpr const double PI_2 = PI*PI;
	constexpr const double PI_3 = PI*PI*PI;

	/** @name Values of the \f$\zeta\f$-function.
	 */
	///@{
	/** \f$\zeta(2)\f$ */
	const double Zeta2 = PI*PI/6.0;
	/** \f$\zeta(2)\f$ */
	const double Zeta3 = 1.2020569031595942854;
	///@}


	/** @name Evolution constants
	 * @todo make these user-configurable
	 */
	///@{
	/** @brief Number of distributions that evolve */
	const uint DISTS = 37;

	/** @brief Number of gauss-legendre abscissae/weights */
	const uint GAUSS_POINTS = 30;

	/** @brief Number of interpolation points to use */
	const uint INTERP_POINTS = 5;

	/** @brief Number of iterations for non-singlet */
	const uint DEFAULT_ITERATIONS = 10;

	/** @brief number of additional iterations for truncated ansatz */
	const uint DEFAULT_TRUNC_IDX = 5;
	///@}



	/** @name Multi-dimensional arrays
	 *  @brief Nicer syntax for large numbers of nested `std::vectors`.
	 */
	///@{
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

	typedef MultiDimVector<double,2>::type Vector2; //!< 2-dim vector
	typedef MultiDimVector<double,3>::type Vector3; //!< 3-dim vector
	typedef MultiDimVector<double,4>::type Vector4; //!< 4-dim vector
	///@}





	/** @brief concept for iterable template parameters, i.e. has begin()/end() funcs
	 */
	template<typename T>
	concept Iterable = requires(T t) {
		{ std::ranges::begin(t) } -> std::input_or_output_iterator;
		{ std::ranges::end(t) } -> std::sentinel_for<decltype(std::ranges::begin(t))>;
	};


	/** @brief concept for arithmetic arguments
	 */
	template <typename T>
	concept Arithmetic = std::integral<T> or std::floating_point<T>;
};


#endif // __COMMON_HPP
