#include "Candia-v2/FuncArrayGrid.hpp"

namespace Candia2
{
	// ArrayGrid::size_type ArrayGrid::size() const noexcept
	// {
	// 	return _cache.size();
	// }

	void ArrayGrid::zero() noexcept
	{
		for (double& _x : _base)
			_x = 0.0;
	}

	double ArrayGrid::operator[](uint idx) const
	{
		return _base[idx];
	}

	double& ArrayGrid::operator[](uint idx)
	{
		return _base[idx];
	}
}
