#include <iostream>
#include <memory>
#include <vector>
#include <filesystem>
#include <fstream>
#include <numeric>
using namespace std;
namespace fs = filesystem;

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Distribution.hpp"
using namespace Candia2;

int main(int argc, char *argv[]) {
	if (argc != 3)
	{
		cerr << "Must provide an order and number of grid points" << endl;
		return 1;
	}

	const uint order = stoi(argv[1]);
	const uint num_grid_points = stoi(argv[2]);
	
	// define the grid points
	// vector<double> xtab{1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab, num_grid_points);

	// initialize the solver, evolve to 100.0
	// use the Les Houche distribution
	const double Qf = 100.0;
	const uint iterations = 8;
	const uint trunc_idx = 7;
	DGLAPSolver solver(order, grid, Qf, iterations, trunc_idx, std::make_unique<LesHouchesDistribution>());

	MultiDimVector<double, 2>::type F = solver.Evolve();

	fs::path filepath{(order == 3) ? "n3lo.dat" : (order == 2) ? "nnlo.dat" : (order == 1) ? "nlo.dat" : "lo.dat"};
	ofstream outfile(filepath.make_preferred());

	outfile << "# using n=" << iterations << " iterations, "
			<< num_grid_points << " grid points, "
			<< "and a truncation index of " << trunc_idx << '\n';
	
	outfile << scientific;

    vector<double> ntab(grid.Size()-1);
	iota(ntab.begin(), ntab.end(), 0);
	array<int, 13> dists{};
	iota(dists.begin(), dists.end(), 0);

	outfile << setw(7) << "x" << ' ';
	vector<uint> ids{26};
	for (const uint j : ids)
	{
		outfile << setw(15) << j << ' ';
	}
	outfile << '\n';
	
	for (uint k=0; k<grid.Size(); k++)
	{
		outfile << setw(7) << setprecision(1) << grid.At(k) << ' ';
		
		outfile << setprecision(8);	
		for (const uint j : ids)
		{
			outfile << setw(15) << F[j][k] << ' ';
		}
		outfile << '\n';
	}
			
	outfile.close();

}
