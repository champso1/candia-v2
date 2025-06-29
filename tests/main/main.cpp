#include <memory>
#include <vector>
using namespace std;

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Distribution.hpp"
using namespace Candia2;

int main() {
	// define the grid points
	vector<double> xtab{1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab, 801);

	// initialize the solver, evolve to 100.0
	// use the Les Houche distribution
	DGLAPSolver solver(0, grid, 100.0, std::make_unique<LesHouchesDistribution>());

	solver.Evolve();

	// solver.OutputDataFileOld("out-LO.dat");
}
