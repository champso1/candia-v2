// lhapdf_roundtrip.cpp
//
// Writes an NNLO LHAPDF grid across the bottom threshold with LHAPDFGrid, reads it back through
// LHAPDF, and checks it against independent evolutions:
//  - every flavour at a Q knot away from the threshold (catches PDG-id mix-ups, e.g. u and d swapped)
//  - just below and just above m_b, where the grid must keep the discontinuity from the matching
//    (b = 0 below; b > 0 and a different gluon above)
//  - alpha_s at a knot
//
// usage: regression_lhapdf_roundtrip <path to infofile.in>

#include "Candia-v2/LHAPDFGrid.hpp"

#include <LHAPDF/LHAPDF.h>

#include <cmath>
#include <cstdio>
#include <memory>

using namespace Candia2;

namespace
{
	constexpr uint order = 2;
	constexpr uint iterations = 6;
	constexpr uint trunc_idx = 4;
	const std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};

	uint failures = 0;

	void check(bool condition, std::string const& what)
	{
		if (!condition) {
			++failures;
			std::fprintf(stderr, "  FAIL: %s\n", what.c_str());
		}
	}

	/** evolves the toy PDFs to Q with nff flavours (nff = 0: chosen by the distribution) */
	std::vector<ArrayGrid> evolveTo(double Q, uint nff, double* alpha_s = nullptr)
	{
		Grid grid(xtab);
		LesHouchesDistribution dist(Q);
		AlphaS alphas(order, dist.Q0(), Q, dist.alpha0(), 1.0);
		alphas.setVFNS(dist.masses(), dist.nfi(), nff == 0 ? dist.nff() : nff);
		if (alpha_s)
			*alpha_s = alphas.pre(alphas.nff()+1);
		DGLAPSolver solver(order, grid, alphas, Q, iterations, trunc_idx, dist, 1.0);
		return solver.evolve();
	}
}

int main(int argc, char* argv[])
{
	if (argc != 2) {
		std::fprintf(stderr, "usage: %s <infofile.in>\n", argv[0]);
		return 2;
	}
	getLogOptions().verbosity = LOG_WARNING;
	LHAPDF::setVerbosity(0);

	constexpr double Qknot = 10.0;
	constexpr double mb = 4.5;

	// write the grid; m_b lies inside the range, so it is split into two subgrids
	{
		LesHouchesDistribution dist(Qknot);
		Grid grid(xtab);
		LHAPDFGrid lhapdf_grid("candia_roundtrip", argv[1], dist, grid, order, iterations, trunc_idx, 1.0);
		lhapdf_grid.evolve({dist.Q0(), 2.0, 3.0, 4.0, 5.0, 6.5, 8.0, Qknot}, {});
		lhapdf_grid.write();
	}
	LHAPDF::pathsPrepend(".");
	std::unique_ptr<LHAPDF::PDF> pdf(LHAPDF::mkPDF("candia_roundtrip", 0));

	Grid grid(xtab);
	auto rel = [](double a, double b) { return std::abs(a - b)/std::max(std::abs(a), std::abs(b)); };

	// 1) every flavour at a knot away from the threshold, plus alpha_s
	{
		double alpha_s{};
		auto F = evolveTo(Qknot, 0, &alpha_s);
		// PDG id -> index in the solver output (0 = g, 1-6 = u d s c b t, 7-12 = antiquarks)
		constexpr std::pair<int, uint> pid_to_index[] = {
			{21, 0}, {2, 1}, {1, 2}, {3, 3}, {4, 4}, {5, 5},
			{-2, 7}, {-1, 8}, {-3, 9}, {-4, 10}, {-5, 11}};
		uint compared = 0;
		for (auto [pid, idx] : pid_to_index) {
			for (uint k : grid.ntab()) {
				const double expected = F[idx][k];
				if (std::abs(expected) < 1e-8)
					continue; // too small for a meaningful relative comparison
				const double got = pdf->xfxQ(pid, grid[k], Qknot);
				++compared;
				check(rel(got, expected) < 1e-6, std::format("pid {} at x={:.0e}, Q={}: LHAPDF={:.10e} solver={:.10e}", pid, grid[k], Qknot, got, expected));
			}
		}
		check(compared > 50, std::format("only {} values compared at Q={}", compared, Qknot));
		check(rel(pdf->alphasQ(Qknot), alpha_s) < 1e-8, std::format("alpha_s(Q={}): LHAPDF={:.10e} solver={:.10e}", Qknot, pdf->alphasQ(Qknot), alpha_s));
	}

	// 2) both sides of the bottom threshold
	{
		auto below = evolveTo(mb, 4); // pre-matching
		auto above = evolveTo(mb, 5); // post-matching
		const double eps = 1e-9;
		for (uint k : grid.ntab()) {
			const double x = grid[k];
			if (x > 0.5)
				continue;
			const double b_below = pdf->xfxQ(5, x, mb*(1.0-eps));
			const double b_above = pdf->xfxQ(5, x, mb*(1.0+eps));
			const double g_below = pdf->xfxQ(21, x, mb*(1.0-eps));
			const double g_above = pdf->xfxQ(21, x, mb*(1.0+eps));
			check(std::abs(b_below) < 1e-12, std::format("b just below m_b at x={:.0e} should be 0, got {:.3e}", x, b_below));
			check(rel(b_above, above[5][k]) < 1e-6, std::format("b just above m_b at x={:.0e}: LHAPDF={:.10e} solver={:.10e}", x, b_above, above[5][k]));
			check(rel(g_below, below[0][k]) < 1e-6, std::format("g just below m_b at x={:.0e}: LHAPDF={:.10e} solver={:.10e}", x, g_below, below[0][k]));
			check(rel(g_above, above[0][k]) < 1e-6, std::format("g just above m_b at x={:.0e}: LHAPDF={:.10e} solver={:.10e}", x, g_above, above[0][k]));
			// at NNLO the matched b is nonzero (and negative at small x)
			check(std::abs(b_above) > 1e-8 && g_above != g_below, std::format("no discontinuity at m_b for x={:.0e}", x));
		}
	}

	std::printf("LHAPDF round trip: %s\n", failures == 0 ? "all checks passed" : "FAILED");
	return failures == 0 ? 0 : 1;
}
