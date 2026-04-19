#include "Candia-v2/Grid.hpp"
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <vector>
#include <fstream>
#include <numeric>
#include <cstdlib>
#include <chrono>
#include <filesystem>
#include <ranges>
using namespace std;
namespace fs = filesystem;

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Distribution.hpp"
using namespace Candia2;
using out_type = std::vector<ArrayGrid>;

static void usage()
{
	cout << "[ERROR] evolve.cpp: Invalid arguments.\n";
	cout << "Usage:\n";
	cout << "-------------------------------------------------------\n";
	cout << "./evolve(.exe) <order> <iterations> <trunc_idx> <kr> [title]\n";
	cout << "    <order>: perturbative order to perform the calculation.\n";
	cout << "    <iterations>: number of total iterations to perform.\n";
	cout << "    <trunc_idx>: number of truncation iterations to perform (for each main iteration!)\n";
	cout << "    <kr>: ratio of mu_R / mu_F.\n";
	cout << "    [title]: optional -- gives a title for the resulting datafile and logfile\n";
	cout << "-------------------------------------------------------\n\n";
	exit(EXIT_FAILURE);
}

static constexpr char const* DATAFILEDIR = "data";

static void outputData(
	out_type const& F, Grid::grid_type const& xtab, Grid const& grid,
	uint order, uint num_grid_points, uint iterations, uint trunc_idx, double kr,
	std::string filename="")
{
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

	// print them out
	for (uint k=0; k<grid.size(); k++){
		outfile << setw(15) << setprecision(8) << grid.at(k) << ' ';
		outfile << setprecision(std::numeric_limits<double>::max_digits10);	
		for (uint j=0; j<DISTS; ++j)
			outfile << F[j][k] << ' ';
		outfile << '\n';
	}
}

int main(int argc, char *argv[]) {
	if (argc != 5 && argc != 6)
		usage();

	const uint order = stoi(argv[1]);
	const uint iterations = stoi(argv[2]);
	const uint trunc_idx = stoi(argv[3]);
	const double kr = stold(argv[4]);
	const double Qf = 100.0;

	std::string datafile_name{};
	if (argc == 6)
		datafile_name = argv[5];

	ostringstream logfile_ss{};
	logfile_ss << ((order == 3) ? "n3lo" : (order == 2) ? "nnlo" : (order == 1) ? "nlo" : "lo");
	logfile_ss << "-i" << iterations << "-t" << trunc_idx << "-r" << setprecision(2) << kr << ".log";
	if (!fs::exists("log")) {
		if (!fs::create_directory("log"))
			log(LOG_ERROR, "evolve.cpp", "Failed to create log output directory");
	}
	fs::path log_path = fs::current_path()/"log"/(datafile_name.empty() ? logfile_ss.str() : datafile_name);
	log_path.replace_extension(".log");
	std::ofstream log_output_file(log_path);

	auto& log_options = getLogOptions();
	log_options.show_debug_messages = false;
	log_options.show_thread_output = true;
	log_options.use_log_output_stream = true;
	log_options.log_output_stream = log_output_file;
	
	vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab, make_grid_filler<GridFillerLogLinQuad>(), {.split_interval = true});
	auto& grid_options = grid.getOptions();
	grid_options.use_alt_mapping = true;
	grid_options.use_gsl_conv_routine = false;
	grid_options.use_gsl_interp_routine = true;
	
	LesHouchesDistribution dist(Qf);
	AlphaS alphas(order, dist.Q0(), Qf, dist.alpha0(), kr);
	auto& alphas_options = alphas.getOptions();
	// alphas_options.use_broken_log_value = true;
	alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
	// alphas.setFFNS(4);

	DGLAPSolver solver(order, grid, alphas, Qf, iterations, trunc_idx, dist, kr);
	auto& dglap_options = solver.getOptions();
	dglap_options.use_truncated_nonsinglet_sol = false;
	dglap_options.use_nnlo_matching_conditions_at_n3lo = true;
	dglap_options.use_n3lo_heavyquark_asymmetry = true;
	dglap_options.use_fortran_n3lo_splitfuncs = false;

	auto t0 = chrono::high_resolution_clock::now();
	auto F = solver.evolve();
	auto tf = chrono::high_resolution_clock::now();
	chrono::duration<double, ratio<60>> mins = tf-t0;
	log(LOG_INFO, "evolve.cpp", "Evolution took {}.", mins);

	datafile_name += ".dat";
	outputData(F, xtab, grid, order, grid.size(), iterations, trunc_idx, kr, datafile_name);

	if (grid.getOptions().use_gsl_conv_routine) {
		auto const& gsl_conv_errors = solver.getGrid().getGSLConvolutionErrors();
		fs::path gsl_conv_errors_log_path("gsl-conv-errors.dat");
		std::ofstream gsl_conv_errors_log_file(gsl_conv_errors_log_path);
		std::ranges::copy(
			gsl_conv_errors | std::views::transform([](auto&& _t){ auto [x,out,res] = _t; return std::format("{} {} {}", x, out, res); }),
			std::ostream_iterator<std::string>(gsl_conv_errors_log_file, "\n"));
	}
}
