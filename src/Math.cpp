#include "Candia-v2/Math.hpp"
#include "Candia-v2/Common.hpp"

namespace Candia2
{

	double factorial(uint x)
	{
		double res = 1.0;
		while (x > 1.0)
			res *= x--;
		return res;
	}	
};
