#include <chrono>
#include <filesystem>
using namespace std;
namespace fs = filesystem;

#include "Candia-v2/Candia.hpp"
using namespace Candia2;
using out_type = std::vector<ArrayGrid>;

static constexpr char const* DATAFILEDIR = "data";

static void outputData(
	out_type const& F, Grid::grid_type const& xtab, Grid const& grid,
	uint order, uint iterations, uint trunc_idx, double kr,
	std::string filename="")
{
	uint num_grid_points = grid.size();
	
	// open the output file, with a filename descriptive of all the provided inputs
	ostringstream outfile_ss{};
	outfile_ss << ((order == 3) ? "n3lo" : (order == 2) ? "nnlo" : (order == 1) ? "nlo" : "lo");
	outfile_ss << "-g" << num_grid_points << "-i" << iterations << "-t" << trunc_idx << "-r" << setprecision(2) << kr << ".dat";
	string outfile_name = outfile_ss.str();
	fs::path datafiledir_path = fs::current_path()/DATAFILEDIR;
	if (!fs::exists(datafiledir_path)) {
		if (!fs::create_directory(datafiledir_path))
			log(LOG_ERROR, "evolve.cpp", "failed to create output directory for datafiles.");
	}
	
	fs::path datafile_path = datafiledir_path/(filename.empty() ? outfile_name : filename);
	ofstream outfile(datafile_path);

	// print a readable header to list the given inputs
	outfile << "# using n=" << iterations << " iterations, "
			<< num_grid_points << " grid points, "
			<< "a truncation index of " << trunc_idx << ", "
			<< "and a scale ratio mu_R^2/mu_F^2 of " << setprecision(2) << kr << '\n';

	// print also the tabulated x values and their corresponding indices in the grid
	outfile << scientific << setprecision(1);
    vector<int> ntab = grid.ntab();
	for (const double x : xtab)
		outfile << x << ' ';
	outfile << '\n';
	for (const int ix : ntab)
		outfile << ix << ' ';
	outfile << '\n';

	outfile << setprecision(std::numeric_limits<double>::max_digits10);
	for (const double x : grid)
		outfile << x << ' ';
	outfile << '\n';

	// print them out
	for (uint k=0; k<grid.size(); k++){
		outfile << setw(15) << setprecision(8) << grid[k] << ' ';
		outfile << setprecision(std::numeric_limits<double>::max_digits10);	
		for (uint j=0; j<F.size(); ++j)
			outfile << F[j][k] << ' ';
		outfile << '\n';
	}
}

int main() {
	uint order = 0;
	uint iterations = 20;
	uint trunc_idx = 15;
	double mur2_muf2 = 1.0;
	double Qf = 100.0;
	std::string datafile_name("qed.dat");


	auto& log_options = getLogOptions();
	log_options.verbosity = LOG_DEBUG;

	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
	Grid grid(xtab, {}, {.num_gauss_points=100});

	QEDDistribution dist(Qf);
	AlphaS alphas(order, dist.Q0(), Qf, dist.alpha0(), mur2_muf2);
	alphas.setFFNS(4);
	AlphaQED alphaem(order, dist.Q0(), Qf, dist.alphaqed0(), mur2_muf2);

	DGLAPSolver solver(order, grid, alphas, alphaem, Qf, iterations, trunc_idx, dist, mur2_muf2);
	solver.getOptions().try_qed = true;
	
	auto t0 = chrono::high_resolution_clock::now();
	auto F = solver.evolveQED();
	auto tf = chrono::high_resolution_clock::now();
	chrono::duration<double, ratio<1>> secs = tf-t0;
	log(LOG_INFO, "qed.cpp", "Evolution took {}.", secs);

	outputData(F, xtab, grid, order, iterations, trunc_idx, mur2_muf2, datafile_name);
}
