#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/Distribution.hpp"
using namespace Candia2;

#include <filesystem>
#include <format>
#include <string>
#include <string_view>
#include <cstdlib>
#include <memory>
#include <numeric>
using namespace std;

using FinalDist = MultiDimVector<double, 2>::type;

void GeneratePlots(std::string_view filename,
				   std::string_view title,
				   std::vector<int> const& dists,
				   std::vector<int> const& ntab,
				   FinalDist const& F,
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
		outfile << grid.At(k) << '\t';
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
	
	// define the grid points
	vector<double> xtab{1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	// vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab, 401);

	// initialize the solver, evolve to 100.0
	// use the Les Houche distribution
	const double Qf = 200.0;
	const uint iterations = 12;
	const uint trunc_idx = 7;
	DGLAPSolver solver(order, grid, Qf, iterations, trunc_idx, make_shared<LesHouchesDistribution>());
	
	FinalDist F = solver.Evolve();
	
	vector<int> dists{0, 1};

	vector<int> ntab(grid.Size()-1);
	iota(ntab.begin(), ntab.end(), 0);

	stringstream filename_ss{};
	filename_ss << ((order == 0) ? "lo" : (order == 1) ? "nlo" : "nnlo");
	filename_ss << ".dat";
	GeneratePlots(filename_ss.str(), "x \t g \t u", dists, ntab, F, grid);
}
