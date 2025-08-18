#include <iostream>
#include <memory>
#include <vector>
#include <filesystem>
#include <fstream>
#include <numeric>
using namespace std;

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Distribution.hpp"
using namespace Candia2;

int main(int argc, char *argv[]) {
	if (argc != 2)
	{
		cerr << "Must provide an order" << endl;
		return 1;
	}

	const uint order = std::stoi(argv[1]);
	
	// define the grid points
	// vector<double> xtab{1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab, 101);

	// initialize the solver, evolve to 100.0
	// use the Les Houche distribution
	const double Qf = 100.0;
	const uint iterations = 8;
	const uint trunc_idx = 7;
	DGLAPSolver solver(order, grid, Qf, iterations, trunc_idx, std::make_unique<LesHouchesDistribution>());

	MultiDimVector<double, 2>::type F = solver.Evolve();

	filesystem::path filepath{"n3lo.dat"};
	ofstream outfile(filepath.make_preferred());
	outfile << scientific << setprecision(9) << setw(15);

    vector<double> ntab(grid.Size()-1);
	iota(ntab.begin(), ntab.end(), 0);
	array<int, 13> dists{};
	iota(dists.begin(), dists.end(), 0);
	
	for (int k : grid.Ntab())
	{
		outfile << grid.At(k) << '\t';

		for (uint j : vector{26, 32, 25})
		{
			outfile << F[j][k] << '\t';
		}
		outfile << '\n';
	}
			
	outfile.close();

}
