#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
using namespace Candia2;
using out_type = std::vector<ArrayGrid>;

#include <chrono>
using namespace std;
namespace fs = filesystem;

static void usage()
{
	cout << "[ERROR] evolve-qed.cpp: Invalid arguments.\n";
	cout << "Usage:\n";
	cout << "-------------------------------------------------------\n";
	cout << "./evolve-qed <iterations> [title]\n";
	cout << "    <iterations>: number of total iterations to perform.\n";
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

	auto num_dists = static_cast<uint>(QEDPartonIndices::COUNT);
	
	// print them out
	for (uint k=0; k<grid.size(); k++){
		outfile << setw(15) << setprecision(8) << grid[k] << ' ';
		outfile << setprecision(std::numeric_limits<double>::max_digits10);	
		for (uint j=0; j<num_dists; ++j)
			outfile << F[j][k] << ' ';
		outfile << '\n';
	}
}

int main(int argc, char *argv[]) {
	if (argc != 2 && argc != 3)
		usage();

	uint order = 0;
	uint iterations = stoi(argv[1]);
	uint trunc_idx = 0;
	double mur2_muf2 = 1.0;
	double Qf = 100.0;

	std::string datafile_name{};
	
	if (argc == 3) {
		datafile_name = argv[2];
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
	log_options.verbosity = LOG_DEBUG;
	log_options.use_log_output_stream = true;
	log_options.log_output_stream = log_output_file;

	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab);

	LesHouchesQED dist(Qf);
	AlphaS alphas(order, dist.Q0(), dist.Qf(), dist.alpha0(), mur2_muf2);
	alphas.setFFNS(4);
	AlphaQED alphaqed(order, dist.Q0(), dist.Qf(), dist.alphaqed0(), mur2_muf2);

	DGLAPSolver solver(order, grid, alphas, alphaqed, Qf, iterations, trunc_idx, dist, mur2_muf2);
	solver.getOptions().try_qed = true;
	
	auto t0 = chrono::high_resolution_clock::now();
	auto F = solver.evolve();
	auto tf = chrono::high_resolution_clock::now();
	chrono::duration<double, ratio<1>> secs = tf-t0;
	log(LOG_INFO, "evolve.cpp", "Evolution took {}.", secs);

	outputData(F, xtab, grid, order, iterations, trunc_idx, mur2_muf2, datafile_name);
}
