#include "Candia-v2/ArrayGrid.hpp"

namespace Candia2
{
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
