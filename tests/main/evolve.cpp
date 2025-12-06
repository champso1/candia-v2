#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>
#include <fstream>
#include <numeric>
#include <cstdlib>
using namespace std;

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Distribution.hpp"
using namespace Candia2;

static void usage()
{
	cout << "[ERROR] Invalid arguments.\n";
	cout << "Usage:\n";
	cout << "-------------------------------------------------------\n";
	cout << "./evolve(.exe) <order> <num_grid_points> <iterations> <trunc_idx> <kr>\n";
	cout << "    <order>: perturbative order to perform the calculation.\n ";
	cout << "    <num_grid_points>: number of grid points to use.\n";
	cout << "    <iterations>: number of total iterations to perform.\n";
	cout << "    <trunc_idx>: number of truncation iterations to perform (for each main iteration!)\n";
	cout << "    <kr>: ratio of mu_R / mu_F.\n";
	cout << "-------------------------------------------------------\n\n";
}

int main(int argc, char *argv[]) {
	if (argc != 6)
	{
		usage();
		exit(EXIT_FAILURE);
	}

	const uint order = stoi(argv[1]);
	const uint num_grid_points = stoi(argv[2]);
	const uint iterations = stoi(argv[3]);
	const uint trunc_idx = stoi(argv[4]);
	const double kr = stold(argv[5]);
	
	// define the tabulated grid points
	vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab, num_grid_points);

	// initialize the solver, evolve to 100.0 GeV
	// use the Les Houche distribution
	const double Qf = 100.0;
	DGLAPSolver solver(order, grid, Qf, iterations, trunc_idx, make_unique<LesHouchesDistribution>(), kr, true);

	// grab the resultant evolved distributions
	auto F = solver.evolve();

	// open the output file, with a filename descriptive of all the provided inputs
	ostringstream outfile_ss{};
	outfile_ss << ((order == 3) ? "n3lo" : (order == 2) ? "nnlo" : (order == 1) ? "nlo" : "lo");
	outfile_ss << "-g" << num_grid_points << "-i" << iterations << "-t" << trunc_idx << "-r" << setprecision(2) << kr << ".dat";
	string outfile_name = outfile_ss.str();
	ofstream outfile(outfile_name);

	// print a readable header to list the given inputs
	outfile << "# using n=" << iterations << " iterations, "
			<< num_grid_points << " grid points, "
			<< "a truncation index of " << trunc_idx << ", "
			<< "and a scale ratio mu_R/mu_F of " << setprecision(2) << kr << '\n';

	// print also the tabulated x values and their corresponding indices in the grid
	outfile << scientific << setprecision(1);
    vector<int> ntab = grid.ntab();
	for (const double x : xtab)
		outfile << x << ' ';
	outfile << '\n';
	for (const int ix : ntab)
		outfile << ix << ' ';
	outfile << '\n';

	// construct all the distributions
	array<int, 13> dists{};
	iota(dists.begin(), dists.end(), 0);

	// print them out
	for (uint k=0; k<grid.size(); k++)
	{
		outfile << setw(15) << setprecision(8) << grid.at(k) << ' ';
		
		outfile << setprecision(8);	
		for (const uint j : dists)
		{
			outfile << setw(15) << F[j][k] << ' ';
		}
		outfile << '\n';
	}
}
