/**
 *  @file ArrayGrid.hpp
 *  @brief Contains the @a ArrayGrid class which is nothing more than an array with some convenience methods
 */

#ifndef __ARRAYGRID_HPP
#define __ARRAYGRID_HPP

#include <concepts>
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
		explicit ArrayGrid(ArrayGrid&& other) = delete; //!< no move constructor

		/**
		 *  @brief Initializes using list-initialization
		 *  @param l an initializer list of values to initialize the array with
		 */
		ArrayGrid(std::initializer_list<value_type> l) : _base(l) {}

		/** Copy assignment operator performs like the copy constructor */
		inline void operator=(ArrayGrid const& other) { _base = other.base(); }
		void operator=(ArrayGrid&& other) = delete; //!< no move assignment operator


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

	
}

#endif // __ARRAYGRID_HPP
