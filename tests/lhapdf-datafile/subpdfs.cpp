#include "Candia-v2/Candia.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"
using namespace Candia2;

#include "LHAPDF/LHAPDF.h"

#include <vector>
#include <fstream>
#include <format>
#include <string>

int main(int argc, char** argv)
{
	LHAPDF::PDF* pdf = LHAPDF::mkPDF("testpdf", 0);

	double Qf = 10.0;
	double Qf2 = Qf*Qf;

	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab);
	
	LesHouchesDistribution dist(Qf);
	double mc = dist.masses(DIST_C);
	AlphaS alphas(3, dist.Q0(), 10.0, dist.alpha0(), 1.0);
	alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
	alphas.calculateThresholdValues();
	double as = alphas.pre(dist.nff()+1);
	double as2 = as*as;
	double as3 = as2*as;

	log(LOG_INFO, "subpdfs.cpp", "using as(nf={}) = {}", dist.nff(), as);

	auto zero_func = [](double,double){ return 0.0; };
	auto a1qg_reg_func = [as](double lm, double nf, double x) {
		auto trunced = ome::AQg_reg.truncate(1);
		return trunced(as, lm, nf, x); };
	auto a2hq_reg_func = [as](double lm, double nf, double x) {
		auto trunced = ome::AQqPS_reg.truncate(2);
		return trunced(as, lm, nf, x); };
	auto a2hg_reg_func = [as](double lm, double nf, double x) {
		auto trunced = ome::AQg_reg.truncate(2);
		return trunced(as, lm, nf, x); };
	auto a3hq_reg_func = [as](double lm, double nf, double x) {
		auto trunced = ome::AQqPS_reg.truncate(3);
		return trunced(as, lm, nf, x); };
	auto a3hg_reg_func = [as](double lm, double nf, double x) {
		auto trunced = ome::AQg_reg.truncate(3);
		return trunced(as, lm, nf, x); };
    
	OpMatElemCustom a1hg(a1qg_reg_func, zero_func, zero_func);
	OpMatElemCustom a2hq(a2hq_reg_func, zero_func, zero_func);
	OpMatElemCustom a2hg(a2hg_reg_func, zero_func, zero_func);
	OpMatElemCustom a3hq(a3hq_reg_func, zero_func, zero_func);
	OpMatElemCustom a3hg(a3hg_reg_func, zero_func, zero_func);
	

	ArrayGrid c(grid.size()), sigma(grid.size()), g(grid.size());
	for (auto&& [i, x] : grid.enumerate()) {
		g[i] = pdf->xfxQ2(21, x, Qf2);
		c[i] = pdf->xfxQ2(4, x, Qf2);
		sigma[i] = 0.0;
		for (uint j=1; j<=4; ++j)
			sigma[i] += pdf->xfxQ2(j, x, Qf2) + pdf->xfxQ2(-j, x, Qf2);
	}

	double L = std::log(std::pow(Qf/mc, 2.0));
	OpMatElemN3LO::update(-L, 4);

	std::ofstream outfile("out.dat");
	outfile << std::scientific << std::setprecision(6);
	for (auto&& [k, x] : grid.enumerate()) {
		double ftilde1 = as*grid.convolution(g, a1hg, k);
		double ftilde2 = as2*(
			grid.convolution(sigma, a2hq, k) +
			grid.convolution(g, a2hg, k)
		);
		double ftilde3 = as3*(
			grid.convolution(sigma, a3hq, k) +
			grid.convolution(g, a3hg, k)
		);
		double ftildennlo = ftilde1 + ftilde2;
		double ftilden3lo = ftilde1 + ftilde2 + ftilde3;
		// double deltaf1 = std::abs(c[k] - ftilde1);
		// double deltaf2 = std::abs(c[k] - ftilde2);
		// double deltaf3 = std::abs(c[k] - ftilde3);
		double deltafnnlo = std::abs(c[k] - ftildennlo);
		double deltafn3lo = std::abs(c[k] - ftilden3lo);

		outfile << x << ' ';
		outfile << c[k] << ' ';
		outfile <<
			ftilde1 << ' ' <<
			ftilde2 << ' ' <<
			ftilde3 << ' ' <<
			ftildennlo << ' ' <<
			ftilden3lo << ' ' <<
			deltafnnlo << ' ' <<
			deltafn3lo << '\n';
			
	}

	delete pdf;
}
