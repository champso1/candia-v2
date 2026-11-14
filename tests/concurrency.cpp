// concurrency.cpp
//
// Runs two different evolutions (which pass through different n_f values) at the same time and
// checks that each gives exactly the same result as when it runs alone. This fails if solvers
// share mutable state, e.g. the n_f of the splitting functions.

#include "Candia-v2/Candia.hpp"

#include <cstdio>
#include <thread>
#include <vector>

using namespace Candia2;

namespace
{
	std::vector<ArrayGrid> run(uint order, double Qf)
	{
		Grid grid({1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9});
		LesHouchesDistribution dist(Qf);
		AlphaS alphas(order, dist.Q0(), Qf, dist.alpha0(), 1.0);
		alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
		DGLAPSolver solver(order, grid, alphas, Qf, 6, 4, dist, 1.0);
		return solver.evolve();
	}

	bool identical(std::vector<ArrayGrid> const& a, std::vector<ArrayGrid> const& b, const char* name)
	{
		if (a.size() != b.size()) {
			std::fprintf(stderr, "%s: different number of distributions\n", name);
			return false;
		}
		for (uint j=0; j<a.size(); ++j) {
			if (a[j].size() != b[j].size()) {
				std::fprintf(stderr, "%s: different grid sizes for distribution %u\n", name, j);
				return false;
			}
			for (uint k=0; k<a[j].size(); ++k) {
				if (a[j][k] != b[j][k]) {
					std::fprintf(stderr, "%s: distribution %u differs at k=%u (%.17g vs %.17g)\n", name, j, k, a[j][k], b[j][k]);
					return false;
				}
			}
		}
		return true;
	}
}

int main()
{
	getLogOptions().verbosity = LOG_WARNING;

	// N3LO to 100 GeV passes through nf=3..5, NNLO to 3 GeV stays at nf=3..4
	auto n3lo_alone = run(3, 100.0);
	auto nnlo_alone = run(2, 3.0);

	std::vector<ArrayGrid> n3lo_together, nnlo_together;
	std::thread t1([&]{ n3lo_together = run(3, 100.0); });
	std::thread t2([&]{ nnlo_together = run(2, 3.0); });
	t1.join();
	t2.join();

	bool ok = identical(n3lo_alone, n3lo_together, "N3LO to 100 GeV")
		& identical(nnlo_alone, nnlo_together, "NNLO to 3 GeV");
	std::printf("concurrent solvers %s their sequential results\n", ok ? "match" : "do NOT match");
	return ok ? 0 : 1;
}
