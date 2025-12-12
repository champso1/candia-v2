#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/FuncArrayGrid.hpp"
using namespace Candia2;

#include <filesystem>
#include <format>
#include <string>
#include <string_view>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <fstream>
#include <iostream>
using namespace std;

void generatePlots(
	std::string_view filename,
	std::string_view title,
	std::vector<int> const& dists,
	std::vector<int> const& ntab,
	std::vector<ArrayGrid> & F,
	Grid const& grid)
{
	string temp{};
	filesystem::path filepath(filename);
		
	ofstream outfile(filepath.make_preferred());
	if (!outfile) {
		temp = format("[PLOT] GeneratePlots(): output file '{}' failed to open", filename);
		cerr << temp << endl;
		exit(EXIT_FAILURE);
	}
		
	
	outfile << "# " << title << '\n';
	for (int k : ntab)
	{
		outfile << scientific << setprecision(9);
		outfile << grid.at(k) << '\t';
		outfile << fixed << setprecision(9);

		for (uint j : dists)
		{
			outfile << F[j][k] << '\t';
		}
		outfile << '\n';
	}
			
	outfile.close();
}

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		cerr << "[PLOT] main(): Must provide an order!" << endl;
		exit(EXIT_FAILURE);
	}
    const uint order = stoi(argv[1]);
	const uint num_grid_points = 1000;
	const uint iterations = 12;
	const uint trunc_idx = 15;
	const double Qf = 100.0;
	const double kr = 1.0;
	
	// define the grid points
	vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid{xtab, num_grid_points};

	std::unique_ptr<LesHouchesDistribution> dist = std::make_unique<LesHouchesDistribution>();
	AlphaS alphas(order, dist->Q0(), dist->alpha0(), Qf, kr);
	alphas.setVFNS(dist->masses(), dist->nfi());

	DGLAPSolver solver(order, grid, alphas, Qf, iterations, trunc_idx, *dist, kr, true);
	
	auto F = solver.evolve();
	vector<int> dists{0, 1};

	vector<int> ntab(grid.size()-1);
	iota(ntab.begin(), ntab.end(), 0);

	stringstream filename_ss{};
	filename_ss << ((order == 0) ? "lo" : (order == 1) ? "nlo" : "nnlo");
	filename_ss << ".dat";
	generatePlots(filename_ss.str(), "x \t g \t u", dists, ntab, F, grid);
}
