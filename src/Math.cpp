#include "Candia-v2/Math.hpp"
#include "Candia-v2/Common.hpp"

#include <cmath>
#include <limits>

namespace Candia2
{

	double Factorial(uint x)
	{
		double res = 1.0;
		while (x > 1.0)
			res *= x--;
		return res;
	}

	
	
};
