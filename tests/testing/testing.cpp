#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
using namespace Candia2;

int main()
{
	auto& options = getLogOptions();
	options.show_debug_messages = true;

	LesHouchesDistribution dist{};
	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab, 500, Grid::LOG_LIN);

	AlphaS alphas(1, dist.Q0(), 100, dist.alpha0(), 1.0);
	alphas.setVFNS(dist.masses(), dist.nfi());

	DGLAPSolver solver(1, grid, alphas, 100.0, 8, 10, dist, 1.0);
	solver.evolve();
}

