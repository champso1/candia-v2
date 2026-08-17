#include "Candia-v2/Candia.hpp"
using namespace Candia2;

#include <vector>
#include <fstream>

int main()
{
	getLogOptions().verbosity = LOG_INFO;

	uint order = 0;
	uint iterations = 10;
	uint trunc_idx = 0;
	double Q0 = 1e-2;
	double Qf = 1e2;
	double logq0 = std::log10(Q0);
	double logqf = std::log10(Qf);
	double mur2_muf2 = 1;
	double a0 = ALPHAQED_MTAU;

	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab);

	LesHouchesQED dist(Qf);
	// LHAPDFDistribution dist(make_lhapdf_pdf("CT18NNLO"), 1.295, 100.0);
	AlphaQED alphaqed(order, dist.Q0(), dist.Qf(), dist.alphaqed0(), mur2_muf2);
	AlphaS alphas(order, dist.Q0(), dist.Qf(), dist.alpha0(), mur2_muf2);
	alphas.setFFNS(4);

	DGLAPSolver solver(order, grid, alphas, alphaqed, Qf, iterations, trunc_idx, dist, mur2_muf2);
	solver.getOptions().try_qed = true;
	auto dists = solver.evolve();
}
