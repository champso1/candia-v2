// api.cpp
//
// Quick checks of specific library behaviour that has regressed before.

// including <cassert> first must work: Candia2 used to define a function named 'assert'
#include <cassert>

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp" // included twice on purpose (include guards)

#include <cmath>
#include <cstdio>
#include <type_traits>

using namespace Candia2;

static_assert(std::is_nothrow_move_constructible_v<ArrayGrid>, "ArrayGrid must be cheaply movable");
static_assert(std::is_nothrow_move_assignable_v<ArrayGrid>, "ArrayGrid must be cheaply movable");
static_assert(!std::is_convertible_v<uint, ArrayGrid>, "ArrayGrid(uint) must stay explicit");

namespace
{
	int failures = 0;

	void check(bool condition, const char* what)
	{
		std::printf("[%s] %s\n", condition ? "PASS" : "FAIL", what);
		if (!condition)
			++failures;
	}

	// max relative error when interpolating a smooth function on the grid
	double interpolationError(Grid& grid)
	{
		auto f = [](double x) { return std::pow(x, 0.3)*std::pow(1.0-x, 3.0); };
		ArrayGrid y(grid.size());
		for (uint k=0; k<grid.size(); ++k)
			y[k] = f(grid[k]);
		double worst = 0.0;
		for (double x : {1e-4, 0.05, 0.5, 0.95, 0.995}) {
			double interp = grid.interpolate(y.view(), x);
			worst = std::max(worst, std::abs(interp - f(x))/f(x));
		}
		return worst;
	}
}

int main()
{
	getLogOptions().verbosity = LOG_WARNING;

	// the grid never contains x=1, even if requested
	{
		Grid grid({1e-5, 1e-3, 0.1, 0.5, 0.9, 1.0});
		check(grid.points().back() < 1.0, "Grid drops x=1 from xtab");
		check(grid.ntab().size() == 5, "Grid keeps the other xtab points");
	}

	// interpolation must use each grid's own size: two grids of different sizes in one process
	{
		Grid small({1e-5, 0.1, 0.9}, GridFillerArgs{.log_size=40, .lin_size=20, .quad_size=10});
		Grid large({1e-5, 0.1, 0.9}, GridFillerArgs{.log_size=200, .lin_size=100, .quad_size=50});
		double err_small = interpolationError(small);
		double err_large = interpolationError(large);
		check(err_small < 1e-3, "interpolation on a small grid is accurate");
		check(err_large < 1e-3, "interpolation on a larger grid created afterwards is accurate");
	}

	// errors are reported as exceptions, even when logging is switched off
	{
		getLogOptions().verbosity = LOG_NONE;
		bool threw = false;
		try {
			AlphaS alphas(2, 1.4, 100.0, 0.35, 1.0);
			alphas.masses(4); // no scheme set
		} catch (Error const&) {
			threw = true;
		}
		check(threw, "errors throw Candia2::Error with LOG_NONE");
		getLogOptions().verbosity = LOG_WARNING;
	}

	// a truncation index below the perturbative order is raised instead of overflowing
	{
		constexpr double Qf = 10.0;
		Grid grid({1e-5, 1e-3, 0.1, 0.5, 0.9});
		LesHouchesDistribution dist(Qf);
		AlphaS alphas(2, dist.Q0(), Qf, dist.alpha0(), 1.0);
		alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
		DGLAPSolver solver(2, grid, alphas, Qf, 4, 0, dist, 1.0);
		auto const& F = solver.evolve();
		bool finite = true;
		for (auto const& dist_k : F)
			for (double v : dist_k)
				finite = finite && std::isfinite(v);
		check(finite, "NNLO evolution with trunc_idx=0 gives finite results");
	}

	// converged, the exact and truncated N3LO non-singlet solutions agree even for charm at x = 0.9, which is ~1e-6
	// of the valence there, so it tests the non-singlet solutions to ~1e-9 (with the beta-function roots that used
	// to be hard-coded to ~1e-6, the exact one was off by 1%); the singlet needs trunc_idx ~ 15 to converge there
	{
		constexpr double Qf = 100.0;
		auto charm = [&](bool truncated) {
			Grid grid({1e-5, 1e-3, 0.1, 0.5, 0.9});
			LesHouchesDistribution dist(Qf);
			AlphaS alphas(3, dist.Q0(), Qf, dist.alpha0(), 1.0);
			alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
			DGLAPSolver solver(3, grid, alphas, Qf, 20, 15, dist, 1.0);
			auto const& F = truncated ? solver.evolveTrunc() : solver.evolve();
			const uint k = grid.ntab()[4];
			return std::pair{F[DIST_C][k], F[DIST_C+6][k]};
		};
		auto [c_exact, cbar_exact] = charm(false);
		auto [c_trunc, cbar_trunc] = charm(true);
		check(std::abs(c_exact/c_trunc - 1.0) < 1e-3 && std::abs(cbar_exact/cbar_trunc - 1.0) < 1e-3,
			"exact and truncated N3LO agree for c and cbar at x = 0.9");
	}

	// subtraction PDFs are computed for the heavy quarks that are active at the final energy
	// (0-9 charm, 10-19 bottom), and are zero for the others
	for (double Qf : {3.0, 10.0}) {
		Grid grid({1e-5, 1e-3, 0.1, 0.5, 0.9});
		LesHouchesDistribution dist(Qf);
		AlphaS alphas(2, dist.Q0(), Qf, dist.alpha0(), 1.0);
		alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
		DGLAPSolver solver(2, grid, alphas, Qf, 4, 4, dist, 1.0);
		auto const& F = solver.evolve();
		auto sub = solver.calculateSubtractionPDFs();
		auto maxAbs = [&sub](uint first, uint last) {
			double out = 0.0;
			for (uint i=first; i<last; ++i)
				for (double v : sub[i])
					out = std::max(out, std::abs(v));
			return out;
		};
		const bool bottom_active = alphas.nff() >= DIST_B;
		const uint k = grid.ntab()[1];
		check(maxAbs(0, 10) > 0.0, "charm subtraction PDFs are computed when charm is active");
		check((maxAbs(10, 20) > 0.0) == bottom_active, "bottom subtraction PDFs are computed only when bottom is active");
		check(std::abs(sub[5][k] - (F[DIST_C][k] - sub[0][k])) < 1e-12, "charm subtraction PDF 5 is c - ftilde1");
	}

	return failures == 0 ? 0 : 1;
}
