#include "Candia-v2/Candia.hpp"
#include "Candia-v2/LHAPDFDistribution.hpp"
using namespace Candia2;

int main()
{
	getLogOptions().show_debug_messages = true;

	const double Q0 = 1.35;
	const double Qf = 100.0;
	const uint order = 3;
	const double mur2_muf2 = 1.0;
	
	LesHouchesDistribution dist(Qf);
	LHAPDFDistribution lhapdf_dist(make_lhapdf_pdf("CT18NLO"), Q0, Qf);
	AlphaS alphas(order, lhapdf_dist.Q0(), Qf, lhapdf_dist.alpha0(), mur2_muf2);
	alphas.setVFNS(lhapdf_dist.masses(), lhapdf_dist.nfi(), lhapdf_dist.nff());
	// alphas.setFFNS(4);
}

