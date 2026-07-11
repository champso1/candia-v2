#include <chrono>
#include <filesystem>
#include <charconv>
#include <string>
namespace fs = std::filesystem;

#include "Candia-v2/Candia.hpp"
using namespace Candia2;
using out_type = std::vector<ArrayGrid>;

static void usage();
static void outputData(
	out_type const& F, std::vector<double> const& xtab, Grid const& grid,
	uint order, uint iterations, uint trunc_idx, double mur2_muf2,
	std::string filename="");


int main(int argc, char *argv[]) {
	if (argc != 6)
		usage();

	constexpr double Qf = 100.0;

	/* ===== cli arguments ===== */
	uint order;       // perturbative order
	uint iterations;  // number of recursive iterations to complete.
	                  // this is the 's' index in the exact ansatz, or the 'n' index in the truncated ansatz
	uint trunc_idx;   // number of terms to keep in the truncated ansatz, corresponding to the 'kappa' or 'k' index
	double mur2_muf2; //
	std::string datafile_name{}; // the name of the output file
	
	order = std::stoi(argv[1]);
	iterations = std::stoi(argv[2]);
	trunc_idx = std::stoi(argv[3]);
	mur2_muf2 = std::stold(argv[4]);
	datafile_name = std::format("{}.dat", argv[7]);
	/* ========== */

	// here we configure how the code logs its output
	// we ensure that anything of "INFO" severity or higher is displayed
	// we can set a logfile where all the output is also sent
	std::ofstream log_output_file("out.log");
	auto& log_options = getLogOptions();
	log_options.verbosity = LOG_INFO;
	log_options.use_log_output_stream = true;
	log_options.log_output_stream = log_output_file;

	// xtab is an array of x-values that we want to ensure are on the grid
	// so that we can grab the values of PDFs at these exact values later
	// without interpolation
	// then we construct the grid of x-values
	// this handles interpolations/convolutions
	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
	Grid grid(xtab);

	// LesHouchesDistribution is a derived class of Distribution
	// these classes handle the initial PDF values as well as the alphas values
	// this distribution implements the LHAPDF toy model
	// where mc = sqrt2, mb = 4.5, mt = 175, Q0 = sqrt2, and as(Q0) = 0.35
	// it also accepts the final energy which, along with Q0, determines the initial/final number of flavors n_f
	LesHouchesDistribution dist(Qf);

	// we then initialize the running coupling, which of course accepts some initial data from th distribution
	// it also takes the perturbative order and mu_R^2/mu_F^2 to perform matching at the quark masss thresholds
	AlphaS alphas(order, dist.Q0(), Qf, dist.alpha0(), mur2_muf2);
	// we can either evolve in the VFNS scheme, in which we need the quark pole mass values to perform the matching
	// as well as the number of the initial and final flavors 
	alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
	// alternatively we can just hard code an n_f value in the FFNS scheme
	// alphas.setFFNS(4);

	// lastly we initialize the solver which accepts most of the above variables/objects
	// this will setup the recursive coefficients (e.g. the D^s_{t,m,n} coefficients at N3LO)
	DGLAPSolver solver(order, grid, alphas, Qf, iterations, trunc_idx, dist, mur2_muf2);
	// there is the option to, at N3LO, change the approximation type of the splitting functions
	// the current status is that they are all available in approximate form, with two envelopes available
	// or the average of the two, given by Imod1 and 2 for each envelope or ImodAvg for the average
	// by default, the P3ns splitting functions utilize the exact expressions,
	// so going back to the approximation can also be done here
	/*
	solver.setP3ApproximationTypes({
			std::make_pair(ExprName::P3nsm, P3ApproxType::ImodAvg),
			std::make_pair(ExprName::P3nsp, P3ApproxType::ImodAvg),
			std::make_pair(ExprName::P3nsv, P3ApproxType::ImodAvg)
		});
	*/

	// lastly, we perform the actual evolution of the coefficients and construction of the final PDFs,
	// which is what is returned.
	// we can either evolve in the exact ansatz (.evolve()) or the truncated ansatz (.evolveTrunc())
	auto F = solver.evolve();
	// auto F = solver.evolveTrunc();

	outputData(F, xtab, grid, order, iterations, trunc_idx, mur2_muf2, datafile_name);
}



static void usage()
{
	log(LOG_ERROR_NOQUIT, "evolve_dglap.cpp", "Invalid arguments.");
	log(LOG_INFO, "evolve_dglap.cpp", "Usage:");
	log(LOG_INFO, "evolve_dglap.cpp", "-------------------------------------------------------");
	log(LOG_INFO, "evolve_dglap.cpp", "./evolve_dglap <order> <iterations> <trunc_idx> <mur2_muf2> <title>");
	log(LOG_INFO, "evolve_dglap.cpp", "    <order>: perturbative order to perform the calculation.");
	log(LOG_INFO, "evolve_dglap.cpp", "    <iterations>: number of total iterations to perform.");
	log(LOG_INFO, "evolve_dglap.cpp", "    <trunc_idx>: number of truncation iterations to perform (for each main iteration!)");
	log(LOG_INFO, "evolve_dglap.cpp", "    <mur2_muf2>: ratio of mu_R^2 / mu_F^2.");
	log(LOG_INFO, "evolve_dglap.cpp", "    <title>: title for the resulting datafile and logfile");
	log(LOG_INFO, "evolve_dglap.cpp", "-------------------------------------------------------\n");
	throw std::runtime_error("invalid cli arguments");
}

static void outputData(
	out_type const& F, Grid::grid_type const& xtab, Grid const& grid,
	uint order, uint iterations, uint trunc_idx, double mur2_muf2,
	std::string filename)
{
	uint num_grid_points = grid.size();
    
	fs::path datafile_path = fs::current_path()/fs::path(filename);
	std::ofstream outfile(datafile_path);

	// print a readable header to list the given inputs
	outfile << "# using n=" << iterations << " iterations, "
			<< num_grid_points << " grid points, "
			<< "a truncation index of " << trunc_idx << ", "
			<< "and a scale ratio mu_R/mu_F of " << std::setprecision(2) << mur2_muf2 << '\n';

	// print also the tabulated x values and their corresponding indices in the grid
	outfile << std::scientific << std::setprecision(1);
	std::vector<int> ntab = grid.ntab();
	for (const double x : xtab)
		outfile << x << ' ';
	outfile << '\n';
	for (const int ix : ntab)
		outfile << ix << ' ';
	outfile << '\n';

	outfile << std::setprecision(std::numeric_limits<double>::max_digits10);
	for (const double x : grid)
		outfile << x << ' ';
	outfile << '\n';

	// print them out
	for (uint k=0; k<grid.size(); k++){
		outfile << std::setw(15) << std::setprecision(8) << grid[k] << ' ';
		outfile << std::setprecision(std::numeric_limits<double>::max_digits10);	
		for (uint j=0; j<DISTS; ++j)
			outfile << F[j][k] << ' ';
		outfile << '\n';
	}
}
