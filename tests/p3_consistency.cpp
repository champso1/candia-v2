// p3_consistency.cpp
//
// Checks that the C++ approximate P3 splitting functions agree with the Fortran (mvv_p3) ones
// for every approximation set, several n_f and several x.
// The residual difference (~1e-7) comes from single-precision literals in the Fortran sources.

#include "Candia-v2/SplittingFn.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace Candia2;

namespace
{
	constexpr double tolerance = 1e-6;

	double relDiff(double a, double b)
	{
		if (a == b)
			return 0.0;
		return std::abs(a - b) / std::max(std::abs(a), std::abs(b));
	}

	const char* imodName(P3ApproxType imod)
	{
		switch (imod) {
			case P3ApproxType::Imod1: return "Imod1";
			case P3ApproxType::Imod2: return "Imod2";
			case P3ApproxType::ImodAvg: return "ImodAvg";
			case P3ApproxType::Exact: return "Exact";
		}
		return "?";
	}

	template <typename TCxx, typename TFortran>
	bool compare(const char* name)
	{
		bool ok = true;
		double worst = 0.0;
		for (auto imod : {P3ApproxType::Imod1, P3ApproxType::Imod2, P3ApproxType::ImodAvg}) {
			for (uint nf : {3u, 4u, 5u}) {
				TCxx cxx(imod);
				TFortran fortran(imod);
				cxx.update(nf, 0.0, 0.0);
				fortran.update(nf, 0.0, 0.0);
				cxx.preCalc();
				fortran.preCalc();

				auto check = [&](const char* part, double x, double a, double b) {
					double d = relDiff(a, b);
					worst = std::max(worst, d);
					if (!(d < tolerance)) {
						ok = false;
						std::fprintf(stderr, "  %s %s nf=%u %s(x=%g): C++=%.12e Fortran=%.12e rel.diff=%.3e\n",
							name, imodName(imod), nf, part, x, a, b, d);
					}
				};
				for (double x : {1e-4, 1e-2, 0.1, 0.4, 0.8, 0.97})
					check("regular", x, cxx.calcRegular(x), fortran.calcRegular(x));
				check("plus", 1.0, cxx.plus(), fortran.plus());
				check("delta", 1.0, cxx.delta(), fortran.delta());
			}
		}
		std::printf("%-6s max rel. diff C++ vs Fortran: %.3e\n", name, worst);
		return ok;
	}
}

int main()
{
	getLogOptions().verbosity = LOG_WARNING;

	bool ok = true;
	ok &= compare<P3nsp, mvv_p3::P3nsp>("P3nsp");
	ok &= compare<P3nsm, mvv_p3::P3nsm>("P3nsm");
	ok &= compare<P3nsv, mvv_p3::P3nsv>("P3nsv");
	ok &= compare<P3qq,  mvv_p3::P3qq>("P3qq");
	ok &= compare<P3qg,  mvv_p3::P3qg>("P3qg");
	ok &= compare<P3gq,  mvv_p3::P3gq>("P3gq");
	ok &= compare<P3gg,  mvv_p3::P3gg>("P3gg");
	return ok ? 0 : 1;
}
