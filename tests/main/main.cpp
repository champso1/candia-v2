#include <memory>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

#include "Candia-v2/Candia.hpp"
#include <Candia-v2/Distribution.hpp>
using namespace Candia2;

int main() {
	// define the grid points
	vector<double> xtab{1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab, 501);

	// pass in our definition for the quark masses/energies
	// we set alpha_s(sqrt(2.0)) = 0.35
	array<double,8> masses{ 0.0 };
	masses[3] = masses[4] = sqrt(2.0);
	masses[5] = 4.5;
	masses[6] = 175.0;
	double Q0 = sqrt(2.0);
	double alpha0 = 0.35;

	// create main alpha_s object
	AlphaS alpha_s(1, 3, Q0, alpha0, masses);

	// initialize the solver, evolve to 200.0
	DGLAPSolver solver(grid, alpha_s, 100.0);

	// set the initial conditions to be the Les Houches toy model
	solver.SetInitialConditions(std::make_shared<LesHouchesDistribution>());

	solver.Evolve();

	solver.OutputDataFile("out.dat");
}
