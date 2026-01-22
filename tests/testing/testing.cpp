#include <vector>

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Grid.hpp"
using namespace Candia2;

int main()
{
	
    getDebugFlag() = true;
	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab, 500, 50, 2);
}
