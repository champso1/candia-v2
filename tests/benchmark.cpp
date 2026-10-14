// benchmark.cpp
//
// Evolves the Les Houches toy PDFs to Q = 100 GeV and compares them with a benchmark table
// (hoppet or pegasus output at the standard xtab points).
//
// usage: regression_benchmark <order> <iterations> <trunc_idx> <mur2_muf2> <truncated: 0/1> <table> <tolerance> [<ffns nf>]
// without <ffns nf> (or with 0) the evolution uses the variable flavour number scheme

#include "Candia-v2/Candia.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace Candia2;

int main(int argc, char* argv[])
{
	if (argc != 8 && argc != 9) {
		std::fprintf(stderr, "usage: %s <order> <iterations> <trunc_idx> <mur2_muf2> <truncated> <table> <tolerance> [<ffns nf>]\n", argv[0]);
		return 2;
	}
	const uint order = std::stoi(argv[1]);
	const uint iterations = std::stoi(argv[2]);
	const uint trunc_idx = std::stoi(argv[3]);
	const double mur2_muf2 = std::stod(argv[4]);
	const bool truncated = std::stoi(argv[5]) == 1;
	const std::string table_path = argv[6];
	const double tolerance = std::stod(argv[7]);
	const uint ffns_nf = (argc == 9) ? std::stoi(argv[8]) : 0;

	getLogOptions().verbosity = LOG_WARNING;

	constexpr double Qf = 100.0;
	const std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
	Grid grid(xtab);
	LesHouchesDistribution dist(Qf);
	AlphaS alphas(order, dist.Q0(), Qf, dist.alpha0(), mur2_muf2);
	if (ffns_nf > 0)
		alphas.setFFNS(ffns_nf);
	else
		alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
	DGLAPSolver solver(order, grid, alphas, Qf, iterations, trunc_idx, dist, mur2_muf2);
	// the hoppet N3LO tables use the approximate non-singlet P3 (the central one of arXiv:1707.08315)
	solver.setP3ApproximationTypes({
		{ExprName::P3nsm, P3ApproxType::ImodAvg},
		{ExprName::P3nsp, P3ApproxType::ImodAvg},
		{ExprName::P3nsv, P3ApproxType::ImodAvg}});
	auto const& F = truncated ? solver.evolveTrunc() : solver.evolve();

	// table columns: x, g, u, d, s, c, b, t, ub, db, sb, cb, bb, tb
	// and the corresponding indices in F (-1 = not compared)
	constexpr int column_to_dist[13] = {0, 1, 2, 3, 4, 5, -1, 7, 8, 9, 10, 11, -1};
	constexpr const char* column_names[13] = {"g", "u", "d", "s", "c", "b", "t", "ub", "db", "sb", "cb", "bb", "tb"};

	std::ifstream table(table_path);
	if (!table) {
		std::fprintf(stderr, "could not open table '%s'\n", table_path.c_str());
		return 2;
	}

	uint row = 0, compared = 0, failures = 0;
	double worst = 0.0;
	std::string line;
	while (std::getline(table, line)) {
		if (line.empty() || line[0] == '#')
			continue;
		std::istringstream ss(line);
		double x{};
		double vals[13]{};
		ss >> x;
		for (double& v : vals)
			ss >> v;
		if (!ss || row >= grid.ntab().size()) {
			std::fprintf(stderr, "malformed table '%s' (row %u)\n", table_path.c_str(), row);
			return 2;
		}
		const uint k = grid.ntab()[row++];

		// x=0.9 is skipped: the heavy flavours are tiny there and the relative difference is meaningless
		if (x > 0.8)
			continue;

		for (int c=0; c<13; ++c) {
			if (column_to_dist[c] < 0 || vals[c] == 0.0)
				continue;
			const double candia = F[column_to_dist[c]][k];
			const double diff = std::abs(candia - vals[c]) / (0.5*std::abs(candia + vals[c]));
			worst = std::max(worst, diff);
			++compared;
			if (!(diff < tolerance)) {
				++failures;
				std::fprintf(stderr, "  %-2s at x=%.0e: candia=%.10e table=%.10e rel.diff=%.3e\n",
					column_names[c], x, candia, vals[c], diff);
			}
		}
	}

	std::printf("order=%u iterations=%u mur2_muf2=%g %s %s: compared %u values, max rel. diff %.3e (tolerance %.1e)\n",
		order, iterations, mur2_muf2, truncated ? "truncated" : "exact", ffns_nf > 0 ? "FFNS" : "VFNS", compared, worst, tolerance);
	if (compared == 0) {
		std::fprintf(stderr, "no values were compared\n");
		return 1;
	}
	return failures == 0 ? 0 : 1;
}
