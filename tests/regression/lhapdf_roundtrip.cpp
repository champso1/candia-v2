// lhapdf_roundtrip.cpp
//
// Writes an LHAPDF grid with LHAPDFGrid, reads it back through LHAPDF, and checks every flavour
// against an independent evolution at one of the Q knots. This catches PDG-id mix-ups
// (e.g. u and d swapped) and formatting problems in the written grid.
//
// usage: regression_lhapdf_roundtrip <path to infofile.in>

#include "Candia-v2/LHAPDFGrid.hpp"

#include <LHAPDF/LHAPDF.h>

#include <cmath>
#include <cstdio>
#include <memory>
#include <numbers>

using namespace Candia2;

int main(int argc, char* argv[])
{
	if (argc != 2) {
		std::fprintf(stderr, "usage: %s <infofile.in>\n", argv[0]);
		return 2;
	}
	getLogOptions().verbosity = LOG_WARNING;
	LHAPDF::setVerbosity(0);

	constexpr uint order = 1;
	constexpr uint iterations = 8;
	constexpr uint trunc_idx = 4;
	constexpr double Q = 10.0;
	const std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
	const std::vector<double> qvals{std::numbers::sqrt2, 2.0, 3.0, 5.0, Q};

	// write the grid
	{
		LesHouchesDistribution dist(qvals.back());
		Grid grid(xtab);
		LHAPDFGrid lhapdf_grid("candia_roundtrip", argv[1], dist, grid, order, iterations, trunc_idx, 1.0);
		lhapdf_grid.evolve(qvals, {});
		lhapdf_grid.write();
	}

	// an independent evolution to the last Q knot
	Grid grid(xtab);
	LesHouchesDistribution dist(Q);
	AlphaS alphas(order, dist.Q0(), Q, dist.alpha0(), 1.0);
	alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
	DGLAPSolver solver(order, grid, alphas, Q, iterations, trunc_idx, dist, 1.0);
	auto const& F = solver.evolve();

	LHAPDF::pathsPrepend(".");
	std::unique_ptr<LHAPDF::PDF> pdf(LHAPDF::mkPDF("candia_roundtrip", 0));

	// PDG id -> index in the solver output (0 = g, 1-6 = u d s c b t, 7-12 = antiquarks)
	constexpr std::pair<int, uint> pid_to_index[] = {
		{21, 0}, {2, 1}, {1, 2}, {3, 3}, {4, 4}, {5, 5},
		{-2, 7}, {-1, 8}, {-3, 9}, {-4, 10}, {-5, 11}};

	uint failures = 0, compared = 0;
	double worst = 0.0;
	for (auto [pid, idx] : pid_to_index) {
		for (uint k : grid.ntab()) {
			double x = grid[k];
			double expected = F[idx][k];
			if (std::abs(expected) < 1e-8)
				continue; // too small for a meaningful relative comparison
			double from_lhapdf = pdf->xfxQ(pid, x, Q);
			double diff = std::abs(from_lhapdf - expected)/std::abs(expected);
			worst = std::max(worst, diff);
			++compared;
			if (!(diff < 1e-6)) {
				++failures;
				std::fprintf(stderr, "  pid %3d at x=%.0e: LHAPDF=%.10e solver=%.10e rel.diff=%.3e\n",
					pid, x, from_lhapdf, expected, diff);
			}
		}
	}

	std::printf("compared %u values read back through LHAPDF, max rel. diff %.3e\n", compared, worst);
	return (failures == 0 && compared > 0) ? 0 : 1;
}
