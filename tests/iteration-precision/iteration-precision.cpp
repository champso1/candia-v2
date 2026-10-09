#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
#include <iomanip>
#include <ios>
#include <iterator>
using namespace Candia2;

#include <vector>
#include <ranges>
#include <fstream>
#include <limits>
#include <sstream>
#include <array>

constexpr static uint order = 3;
constexpr static uint trunc_idx = 15;
constexpr static double Qf = 100.0;
constexpr static double mur2_muf2 = 1.0;

constexpr static uint max_iters = 18;

constexpr static std::array<uint,4> ntabs{0, 4, 6, 8};

std::vector<double> run_dglap(uint iterations)
{
	getLogOptions().verbosity = LOG_WARNING;
	std::vector<double> xtab{1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab);

	LesHouchesDistribution dist(Qf);
	AlphaS alphas(order, dist.Q0(), Qf, dist.alpha0(), mur2_muf2);
	alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
	// alphas.setFFNS(4);

	DGLAPSolver solver(order, grid, alphas, Qf, iterations, trunc_idx, dist, mur2_muf2);
	solver.setP3ApproximationTypes({
			std::make_pair(ExprName::P3nsm, P3ApproxType::ImodAvg),
			std::make_pair(ExprName::P3nsp, P3ApproxType::ImodAvg),
			std::make_pair(ExprName::P3nsv, P3ApproxType::ImodAvg)
		});
    
	auto F = solver.evolveTrunc();
	uint x_01_idx = 4;
	getLogOptions().verbosity = LOG_INFO;
	std::vector<double> vals{};
	for (uint ntab : ntabs)
		vals.emplace_back(F[0][grid.ntab()[ntab]]);
	return vals;
}

using data_ret_type = std::pair<std::vector<uint>, std::vector<std::vector<double>>>;

data_ret_type generate_data()
{
	
	auto iters = std::views::iota(uint{1}, max_iters);
	std::vector<std::vector<double>> gluon_vals{};
	for (uint i=0; i<iters.size(); ++i) {
		log(LOG_INFO, "generate_data()", "iter {}/{}", i+1, iters.size());
		gluon_vals.emplace_back(run_dglap(iters[i]));
	}

	std::ofstream outfile("data.dat");
	outfile << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
	for (uint i=0; i<iters.size(); ++i) {
		outfile << iters[i];
		for (uint j=0; j<ntabs.size(); ++j)
			outfile << ' ' << gluon_vals[i][j];
		outfile << '\n';
	}

	return std::make_pair(std::vector<uint>(iters.begin(), iters.end()), gluon_vals);
}

data_ret_type read_data()
{
	std::vector<uint> iters{};
	std::vector<std::vector<double>> all_gluon_vals{};
	std::ifstream infile("data.dat");
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
	auto[iters, gluon_vals] = generate_data();
	// auto[iters, gluon_vals] = read_data();

	for (uint j=0; j<gluon_vals[0].size(); ++j) {
		double prev=gluon_vals[0][j], cur;
		for (uint i=1; i<iters.size(); ++i) {
			cur = gluon_vals[i][j];
			gluon_vals[i][j] = std::abs(cur-prev)/std::abs((cur+prev)/2.0)*100.0;
			prev = cur;
		}
	}

	auto iters_fixed = iters | std::views::drop(1);
	auto gluon_vals_fixed = gluon_vals | std::views::drop(1);
	std::ofstream outfile("out.dat");
	outfile << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
	for (uint i=0; i<iters_fixed.size(); ++i) {
		outfile << iters_fixed[i];
		for (uint j=0; j<ntabs.size(); ++j)
			outfile << ' ' << gluon_vals_fixed[i][j];
		outfile << '\n';
	}
}
