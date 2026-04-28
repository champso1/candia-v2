#include "Candia-v2/Grid.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"
using namespace Candia2;

#include "LHAPDF/LHAPDF.h"
using namespace LHAPDF;

int main()
{
	PDF* pdf = mkPDF("CT18NLO", 0);
	std::vector<int> pids = pdf->flavors();
 
	vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	GridFillerLogLinQuad filler{};
	Grid grid(xtab, filler, {});

	double q = 10.0;
	double q2 = q*q;
	double as = pdf->alphasQ2(q2)/(4.0*PI);
	double as2 = as*as;
	double mb = pdf->quarkMass(5);
	double mb2 = mb*mb;
	
	getLogOptions().show_debug_messages = true;
	log(LOG_DEBUG, "lhapdf-datafile.cpp", "Using:");
	log(LOG_DEBUG, "lhapdf-datafile.cpp", "  q  ={}", q);
	log(LOG_DEBUG, "lhapdf-datafile.cpp", "  q2 ={}", q2);
	log(LOG_DEBUG, "lhapdf-datafile.cpp", "  as ={}", as);
	log(LOG_DEBUG, "lhapdf-datafile.cpp", "  as2={}", as2);
	log(LOG_DEBUG, "lhapdf-datafile.cpp", "  mb ={}", mb);

	double lm = std::log(mb2/q2);
	uint nf = 5;
	OpMatElem::update(lm, nf);

	auto zero_func = [](double,double,double){ return 0.0; };
	auto delta_zero_func = [](double,double){ return 0.0; };

	auto a1hg_reg_func = [as,lm,nf](double, double, double x) -> double {
		auto trunced = ome::AQg_reg.truncate(1);
		return trunced(as, lm, nf, x); };
	OpMatElemCustom a1hg(a1hg_reg_func, zero_func, delta_zero_func);

	auto a2hq_reg_func = [as](double lm, double nf, double x) {
		auto trunced = ome::AQqPS_reg.truncate(2);
		return trunced(as, lm, nf, x); };
	OpMatElemCustom a2hq(a2hq_reg_func, zero_func, delta_zero_func);
		
	auto a2hg_reg_func = [as](double lm, double nf, double x) {
		auto trunced = ome::AQg_reg.truncate(2);
		return trunced(as, lm, nf, x); };
	OpMatElemCustom a2hg(a2hg_reg_func, zero_func, delta_zero_func);
	
	ArrayGrid b(grid.size()), g(grid.size()), sigma(grid.size());
	for (uint k=0; k<grid.size(); ++k) {
		b[k] = pdf->xfxQ2(5, grid[k], q2);
		g[k] = pdf->xfxQ2(21, grid[k], q2);
		sigma[k] = 0.0;
		for (int nf=1; nf<=5; ++nf)
			sigma[k] += pdf->xfxQ2(nf, grid[k], q2) + pdf->xfxQ2(-nf, grid[k], q2);
	}

	std::cout << "before printing\n";

	std::ofstream outfile("out.dat");
	ArrayGrid ftilde1(grid.size()), ftilde2(grid.size()), ftildeNLO(grid.size());
	for (uint k=0; k<grid.size(); ++k) {
	    ftilde1[k] = grid.convolution(g, a1hg, k);
		ftilde2[k] = as2*(
			grid.convolution(sigma, a2hq, k) +
			grid.convolution(g, a2hg, k));
		ftildeNLO[k] = ftilde1[k] + ftilde2[k];

		outfile << grid[k] << ' '
				<< b[k] << ' '
				<< ftilde1[k] << ' '
				<< ftilde2[k] << ' '
				<< ftildeNLO[k] << '\n';
	}
 
	delete pdf;
	return 0;
}
