#include "Candia-v2/Expression.hpp"
#include <chrono>
#include <filesystem>
using namespace std;
namespace fs = filesystem;

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/LHAPDFDistribution.hpp"
using namespace Candia2;
using out_type = std::vector<ArrayGrid>;

#include "yaml-cpp/yaml.h"

static void usage()
{
	cout << "[ERROR] evolve.cpp: Invalid arguments.\n";
	cout << "Usage:\n";
	cout << "-------------------------------------------------------\n";
	cout << "./evolve(.exe) <order> <iterations> <trunc_idx> <mur2_muf2> <use_trunc> <debug> [title]\n";
	cout << "    <order>: perturbative order to perform the calculation.\n";
	cout << "    <iterations>: number of total iterations to perform.\n";
	cout << "    <trunc_idx>: number of truncation iterations to perform (for each main iteration!)\n";
	cout << "    <mur2_muf2>: ratio of mu_R^2 / mu_F^2.\n";
	cout << "    <use_trunc>: 0=use exact, 1=use truncated";
	cout << "    <debug>: 0=do not show debug messages, 1=show debug messages";
	cout << "    [title]: optional -- gives a title for the resulting datafile and logfile\n";
	cout << "-------------------------------------------------------\n\n";
	throw std::runtime_error("invalid cli arguments");
}

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

	outfile << setprecision(std::numeric_limits<double>::max_digits10);
	for (const double x : grid)
		outfile << x << ' ';
	outfile << '\n';

	// print them out
	for (uint k=0; k<grid.size(); k++){
		outfile << setw(15) << setprecision(8) << grid[k] << ' ';
		outfile << setprecision(std::numeric_limits<double>::max_digits10);	
		for (uint j=0; j<DISTS; ++j)
			outfile << F[j][k] << ' ';
		outfile << '\n';
	}
}

struct CfgFileArgs final
{
	uint order;
	uint iterations;
	uint trunc_idx;
	double mur2_muf2;
	double Qf;
	bool use_trunc;
	bool debug;
	std::string title;
};

static CfgFileArgs read_options_from_yaml()
{
	YAML::Node configfile = YAML::LoadFile("evolve.config.yaml");
    return {
		.order = configfile["order"].as<uint>(),
		.iterations = configfile["iterations"].as<uint>(),
		.trunc_idx = configfile["trunc_idx"].as<uint>(),
		.mur2_muf2 = configfile["mur2_muf2"].as<double>(),
		.Qf = configfile["Qf"].as<double>(),
		.use_trunc = configfile["use_trunc"].as<bool>(),
		.debug = configfile["debug"].as<bool>(),
		.title = configfile["title"].as<std::string>(),
	};
}

int main(int argc, char *argv[]) {
	if (argc != 7 && argc != 8 && argc != 1)
		usage();

	uint order;
	uint iterations;
	uint trunc_idx;
	double mur2_muf2;
	bool use_trunc;
	bool debug;
	double Qf;

	std::string datafile_name{};

	if (argc == 1) {
		auto args = read_options_from_yaml();
		order = args.order;
		iterations = args.iterations;
		trunc_idx = args.trunc_idx;
		mur2_muf2 = args.mur2_muf2;
		use_trunc = args.use_trunc;
		debug = args.debug;
		Qf = args.Qf;
		datafile_name = args.title + ".dat";
	} else {
		order = stoi(argv[1]);
		iterations = stoi(argv[2]);
		trunc_idx = stoi(argv[3]);
		mur2_muf2 = stold(argv[4]);
		use_trunc = stoi(argv[5]) == 1;
		debug = stoi(argv[6]) == 1;
		Qf = 100.0;
	}
	
	if (argc == 8) {
		datafile_name = argv[7];
		datafile_name += ".dat";
	}

	ostringstream fileprefix{};
	fileprefix << ((order == 3) ? "n3lo" : (order == 2) ? "nnlo" : (order == 1) ? "nlo" : "lo");
	fileprefix << "-i" << iterations << "-t" << trunc_idx << "-r" << setprecision(2) << mur2_muf2;

	if (argc == 7)
		datafile_name = fileprefix.str() + ".dat";

	if (!fs::exists("log")) {
		if (!fs::create_directory("log"))
			log(LOG_ERROR, "evolve.cpp", "Failed to create log output directory");
	}
	std::string logfile = fileprefix.str() + ".log";
	fs::path log_path = fs::current_path()/"log"/logfile;
	log_path.replace_extension(".log");
	std::ofstream log_output_file(log_path);

	auto& log_options = getLogOptions();
	log_options.verbosity = debug ? LOG_DEBUG : LOG_INFO;
	log_options.use_log_output_stream = true;
	log_options.log_output_stream = log_output_file;

	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab);

	LesHouchesDistribution dist(Qf);
	// LHAPDFDistribution dist(make_lhapdf_pdf("CT18NNLO"), 1.295, 100.0);
	AlphaS alphas(order, dist.Q0(), Qf, dist.alpha0(), mur2_muf2);
	alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
	// alphas.setFFNS(4);

	DGLAPSolver solver(order, grid, alphas, Qf, iterations, trunc_idx, std::move(dist), mur2_muf2);
	auto& dglap_options = solver.getOptions();
	dglap_options.disable_heavy_flavor_matching = false;
	dglap_options.use_nnlo_matching_conditions_at_n3lo = false;
	dglap_options.use_n3lo_heavyquark_asymmetry = true;
	dglap_options.use_fortran_n3lo_splitfuncs = false;

	auto t0 = chrono::high_resolution_clock::now();
	auto F = use_trunc ? solver.evolveTrunc() : solver.evolve();
	auto tf = chrono::high_resolution_clock::now();
	chrono::duration<double, ratio<1>> secs = tf-t0;
	log(LOG_INFO, "evolve.cpp", "Evolution took {}.", secs);

	outputData(F, xtab, grid, order, iterations, trunc_idx, mur2_muf2, datafile_name);
}
