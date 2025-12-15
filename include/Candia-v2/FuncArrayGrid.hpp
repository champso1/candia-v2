#ifndef __FUNCARRGRID_HPP
#define __FUNCARRGRID_HPP

#include "Candia-v2/Common.hpp"

namespace Candia2
{
	class ArrayGrid final
	{
	private:
		using base_type = std::vector<double>;
		using size_type = std::size_t;

		base_type _base{};
	public:
		explicit ArrayGrid(size_type size)
			: _base(size, 0.0)
		{}
		ArrayGrid(base_type const& points)
			: _base{points}
		{}

		ArrayGrid(ArrayGrid const& other)
			: _base(other._base)
		{}
		inline void operator=(ArrayGrid const& other)
		{
			_base = other._base;
		}

		ArrayGrid(ArrayGrid&& other) = delete;
		void operator=(ArrayGrid&& other) = delete;

		~ArrayGrid() = default;

		inline base_type  const& base() const noexcept { return _base; }
		void zero() noexcept;
		
		double operator[](uint idx) const;  // accessor for points on the grid (const)
		double& operator[](uint idx);       // ditto (non-const)
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
