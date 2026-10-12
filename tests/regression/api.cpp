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

	return failures == 0 ? 0 : 1;
}
