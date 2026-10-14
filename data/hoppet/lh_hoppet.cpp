// lh_hoppet.cpp
//
// Generates the hoppet tables in this directory: the Les Houches toy PDFs evolved from Q0 = sqrt(2) GeV
// to Q = 100 GeV with alpha_s(sqrt(2) GeV) = 0.35, written at the standard x values.
// Not built by CMake; it needs hoppet (written for hoppet 1.2.0's v1 interface):
//
//   g++ -O2 $(hoppet-config --cxxflags) lh_hoppet.cpp -o lh_hoppet $(hoppet-config --libs) -lgfortran
//   ./lh_hoppet <nloop> <mur2_muf2> <ffn nf, or 0 for the VFNS> <output file>
//
// e.g. out-nlo-ffn4-mu05.dat is `./lh_hoppet 2 0.5 4 out-nlo-ffn4-mu05.dat`.
// Note that hoppet's muR_Q is mu_R/mu_F, the square root of mur2_muf2.

#include "hoppet_v1.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>

static void lesHouches(const double& x, const double& Q, double* f)
{
	(void)Q;
	for (int i=0; i<13; ++i)
		f[i] = 0.0;
	const double uv = 5.1072*std::pow(x, 0.8)*std::pow(1.0-x, 3);
	const double dv = 3.06432*std::pow(x, 0.8)*std::pow(1.0-x, 4);
	const double dbar = 0.1939875*std::pow(x, -0.1)*std::pow(1.0-x, 6);
	const double ubar = dbar*(1.0-x);
	// index 6+i holds flavour i in the LHAPDF numbering (0 = gluon)
	f[6] = 1.7*std::pow(x, -0.1)*std::pow(1.0-x, 5);
	f[6+2] = uv + ubar;
	f[6-2] = ubar;
	f[6+1] = dv + dbar;
	f[6-1] = dbar;
	f[6+3] = f[6-3] = 0.2*(ubar + dbar);
}

int main(int argc, char** argv)
{
	if (argc != 5) {
		std::fprintf(stderr, "usage: %s <nloop> <mur2_muf2> <ffn nf, or 0 for the VFNS> <output file>\n", argv[0]);
		return 2;
	}
	const int nloop = std::atoi(argv[1]);
	const double muR_Q = std::sqrt(std::atof(argv[2]));
	const int ffn = std::atoi(argv[3]);

	const double dy = 0.05;
	hoppetStartExtended(12.0, dy, 1.0, 28000.0, dy/4.0, nloop, -6, factscheme_MSbar);
	if (ffn > 0)
		hoppetSetFFN(ffn);
	else
		hoppetSetPoleMassVFN(std::sqrt(2.0), 4.5, 175.0);
	hoppetEvolve(0.35, std::sqrt(2.0), nloop, muR_Q, lesHouches, std::sqrt(2.0));

	FILE* out = std::fopen(argv[4], "w");
	if (!out)
		return 1;
	std::fprintf(out, "#   x        g           u           d           s           c            b          t           ub          db          sb          cb          bb          tb\n");
	// output order: g u d s c b t ub db sb cb bb tb
	constexpr int columns[13] = {6, 8, 7, 9, 10, 11, 12, 4, 5, 3, 2, 1, 0};
	for (double x : {1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9}) {
		double f[13];
		hoppetEval(x, 100.0, f);
		std::fprintf(out, "%.1E", x);
		for (int i : columns)
			std::fprintf(out, " %.10E", f[i]);
		std::fprintf(out, "\n");
	}
	return std::fclose(out) == 0 ? 0 : 1;
}
