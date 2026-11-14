#include "Candia-v2/Candia.hpp"
using namespace Candia2;

#include <vector>
#include <ranges>
#include <fstream>
#include <limits>
#include <sstream>
#include <array>

constexpr static uint order = 3;
constexpr static uint trunc_idx_def = 15;
constexpr static uint iterations_def = 15;
constexpr static double Qf = 100.0;
constexpr static double mur2_muf2 = 1.0;

constexpr static uint max_iters_trunc = 18;

constexpr static std::array<uint,4> ntabs{0, 4, 6, 8};

struct VerbosityHandle final
{
private:
	int _hold;
public:
	VerbosityHandle(){
		_hold = getLogOptions().verbosity;
		getLogOptions().verbosity = LOG_WARNING;
	}
	~VerbosityHandle() { getLogOptions().verbosity = _hold; }
};

std::vector<double> run_dglap_iters(uint iterations)
{
	VerbosityHandle verbosity_handle{};
	
	std::vector<double> xtab{1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
	Grid grid(xtab);

	LesHouchesDistribution dist(Qf);
	AlphaS alphas(order, dist.Q0(), Qf, dist.alpha0(), mur2_muf2);
	alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
	// alphas.setFFNS(4);

	DGLAPSolver solver(order, grid, alphas, Qf, iterations, trunc_idx_def, dist, mur2_muf2);
    
	auto F = solver.evolve();
	std::vector<double> vals{};
	for (uint ntab : ntabs)
		vals.emplace_back(F[0][grid.ntab()[ntab]]);
	return vals;
}

std::vector<double> run_dglap_trunc(uint trunc)
{
	VerbosityHandle verbosity_handle{};
	
	std::vector<double> xtab{1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
	Grid grid(xtab);

	LesHouchesDistribution dist(Qf);
	AlphaS alphas(order, dist.Q0(), Qf, dist.alpha0(), mur2_muf2);
	alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
	// alphas.setFFNS(4);

	DGLAPSolver solver(order, grid, alphas, Qf, iterations_def, trunc, dist, mur2_muf2);
    
	auto F = solver.evolveTrunc();
	std::vector<double> vals{};
	for (uint ntab : ntabs)
		vals.emplace_back(F[0][grid.ntab()[ntab]]);
	return vals;
}

using data_ret_type = std::pair<std::vector<uint>, std::vector<std::vector<double>>>;

std::pair<data_ret_type,data_ret_type> generate_data()
{
	
	auto iters = std::views::iota(uint{1}, max_iters_trunc);
	log(LOG_INFO, "generate_data()", "using iter vals of: {}", vec_to_str2(iters));
	std::vector<std::vector<double>> gluon_vals_iters{};
	for (uint i=0; i<iters.size(); ++i) {
		log(LOG_INFO, "generate_data (iters) ()", "iter {}/{}", i+1, iters.size());
		gluon_vals_iters.emplace_back(run_dglap_iters(iters[i]));
	}

	auto trunc = std::views::iota(uint{3}, max_iters_trunc);
	log(LOG_INFO, "generate_data()", "using trunc vals of: {}", vec_to_str2(trunc));
	std::vector<std::vector<double>> gluon_vals_trunc{};
	for (uint i=0; i<trunc.size(); ++i) {
		log(LOG_INFO, "generate_data (trunc) ()", "iter {}/{}", i+1, trunc.size());
		gluon_vals_trunc.emplace_back(run_dglap_trunc(trunc[i]));
	}

	{
		std::ofstream outfile("data-iters.dat");
		outfile << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
		for (uint i=0; i<iters.size(); ++i) {
			outfile << iters[i];
			for (uint j=0; j<ntabs.size(); ++j)
				outfile << ' ' << gluon_vals_iters[i][j];
			outfile << '\n';
		}
	}
	{
		std::ofstream outfile("data-trunc.dat");
		outfile << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
		for (uint i=0; i<trunc.size(); ++i) {
			outfile << trunc[i];
			for (uint j=0; j<ntabs.size(); ++j)
				outfile << ' ' << gluon_vals_trunc[i][j];
			outfile << '\n';
		}
	}

	std::vector<uint> iters_vec(iters.begin(), iters.end());
	std::vector<uint> trunc_vec(trunc.begin(), trunc.end());

	return std::make_pair(
		std::make_pair(iters_vec, gluon_vals_iters),
		std::make_pair(trunc_vec, gluon_vals_trunc));
}

data_ret_type read_data(std::filesystem::path const& p)
{
	if (p.empty() || !std::filesystem::exists(p))
		log(LOG_ERROR, "iteration-precision.cpp", "filepath {} failed to open", p.string());
	
	std::vector<uint> iters{};
	std::vector<std::vector<double>> all_gluon_vals{};
	std::ifstream infile(p);
	std::string line;
	while (std::getline(infile, line)) {
		std::istringstream ss(line);
		uint iter;
		ss >> iter;
		std::vector<double> gluon_vals(std::istream_iterator<double>(ss), std::istream_iterator<double>{});
		iters.emplace_back(iter);
		all_gluon_vals.emplace_back(gluon_vals);
	}

	return std::make_pair(iters, all_gluon_vals);
}

int main()
{
	auto[iters_vals, trunc_vals] = generate_data();
	auto[iters, gluon_vals_iters] = iters_vals;
	auto[trunc, gluon_vals_trunc] = trunc_vals;
	// auto[iters, gluon_vals_iters] = read_data("data-iters.dat");
	// auto[trunc, gluon_vals_trunc] = read_data("data-trunc.dat");

	for (uint j=0; j<gluon_vals_iters[0].size(); ++j) {
		double prev=gluon_vals_iters[0][j], cur;
		for (uint i=1; i<iters.size(); ++i) {
			cur = gluon_vals_iters[i][j];
			gluon_vals_iters[i][j] = std::abs(cur-prev)/std::abs((cur+prev)/2.0)*100.0;
			prev = cur;
		}
	}
	for (uint j=0; j<gluon_vals_trunc[0].size(); ++j) {
		double prev=gluon_vals_trunc[0][j], cur;
		for (uint i=1; i<trunc.size(); ++i) {
			cur = gluon_vals_trunc[i][j];
			gluon_vals_trunc[i][j] = std::abs(cur-prev)/std::abs((cur+prev)/2.0)*100.0;
			prev = cur;
		}
	}

	{
		auto iters_fixed = iters | std::views::drop(1);
		auto gluon_vals_iters_fixed = gluon_vals_iters | std::views::drop(1);
		std::ofstream outfile("out-iters.dat");
		outfile << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
		for (uint i=0; i<iters_fixed.size(); ++i) {
			outfile << iters_fixed[i];
			for (uint j=0; j<ntabs.size(); ++j)
				outfile << ' ' << gluon_vals_iters_fixed[i][j];
			outfile << '\n';
		}
	}
	{
		auto trunc_fixed = trunc | std::views::drop(1);
		auto gluon_vals_trunc_fixed = gluon_vals_trunc | std::views::drop(1);
		std::ofstream outfile("out-trunc.dat");
		outfile << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
		for (uint i=0; i<trunc_fixed.size(); ++i) {
			outfile << trunc_fixed[i];
			for (uint j=0; j<ntabs.size(); ++j)
				outfile << ' ' << gluon_vals_trunc_fixed[i][j];
			outfile << '\n';
		}
	}
}
