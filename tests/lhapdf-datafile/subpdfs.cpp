#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"
#include <LHAPDF/Config.h>
using namespace Candia2;

#include "LHAPDF/LHAPDF.h"

#include <vector>
#include <fstream>
#include <format>
#include <string>

int main(int argc, char** argv)
{
	getLogOptions().verbosity = LOG_INFO;
	LHAPDF::setVerbosity(0);
	LHAPDF::PDF* pdf = LHAPDF::mkPDF("testpdf", 0);

	double Qf = 4.0;
	double Qf2 = Qf*Qf;

	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(
		xtab,
		{.min=1.0e-5, .log_size=500, .lin_size=200, .quad_size=100, .pivot1=0.1, .pivot2=0.9});
	
	LesHouchesDistribution dist(Qf);
	double mc = dist.masses(DIST_C);
	double mb = dist.masses(DIST_B);
	double mc2 = mc*mc;
	double mb2 = mb*mb;
	AlphaS alphas(3, dist.Q0(), dist.Qf(), dist.alpha0(), 1.0);
	alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
	alphas.calculateThresholdValues();
	double as = alphas.pre(dist.nff()+1)/(4.0*PI);
	// double as = pdf->alphasQ2(Qf2)/(4.0*PI);
	// as = 0.118/(4.0*PI);
	// double mb = pdf->quarkMass(5);
	double as2 = as*as;
	double as3 = as2*as;
	uint nf = 5;

	double L = std::log(mc2/Qf2);
	OpMatElemN3LO::update(L, nf);

	log(LOG_INFO, "subpdfs.cpp", "using alphas(nf={})/4pi = {}, mc={}, L=log(mc2/mu2)={}", nf, as, mc, L);

	auto zero_func = [](double,double){ return 0.0; };
	auto a1qg_reg_func = [as](double lm, double nf, double x) {
		auto omg_reg = ome::AQg_reg[1];
		return omg_reg(lm, nf, x); };
	auto a2hq_reg_func = [as](double lm, double nf, double x) {
		auto omg_reg = ome::AQqPS_reg[2];
		return omg_reg(lm, nf, x); };
	auto a2hg_reg_func = [as](double lm, double nf, double x) {
		auto omg_reg = ome::AQg_reg[2];
		return omg_reg(lm, nf, x); };
	auto a3hq_reg_func = [as](double lm, double nf, double x) {
		auto omg_reg = ome::AQqPS_reg[3];
		return omg_reg(lm, nf, x); };
	auto a3hg_reg_func = [as](double lm, double nf, double x) {
		auto omg_reg = ome::AQg_reg[3];
		return omg_reg(lm, nf, x); };
    
	OpMatElemCustom a1hg(a1qg_reg_func, zero_func, zero_func);
	OpMatElemCustom a2hq(a2hq_reg_func, zero_func, zero_func);
	OpMatElemCustom a2hg(a2hg_reg_func, zero_func, zero_func);
	OpMatElemCustom a3hq(a3hq_reg_func, zero_func, zero_func);
	OpMatElemCustom a3hg(a3hg_reg_func, zero_func, zero_func);
	

	ArrayGrid
		c(grid.size()),
		b(grid.size()),
		sigma(grid.size()),
		g(grid.size());
	for (auto&& [i, x] : grid.enumerate()) {
		c[i] = pdf->xfxQ2(4, x, Qf2);
		b[i] = pdf->xfxQ2(5, x, Qf2);
		g[i] = pdf->xfxQ2(21, x, Qf2);
		sigma[i] = 0.0;
		for (uint j=1; j<=nf; ++j)
			sigma[i] += pdf->xfxQ2(j, x, Qf2) + pdf->xfxQ2(-j, x, Qf2);
	}

	std::ofstream outfile("out-c-nf5.dat");
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
		double ftildenlo  = ftilde1 + ftilde2;
		double ftildennlo = ftilde1 + ftilde2 + ftilde3;
		double deltaf1    = c[k] - ftilde1;
		double deltaf2    = c[k] - ftilde2;
		double deltaf3    = c[k] - ftilde3;
		double deltafnlo  = c[k] - ftildenlo;
		double deltafnnlo = c[k] - ftildennlo;

		outfile
/*1*/			<< x << ' '
/*2*/			<< c[k] << ' ' <<
/*3*/			ftilde1 << ' ' <<
/*4*/			ftilde2 << ' ' <<
/*5*/			ftilde3 << ' ' <<
/*6*/			ftildenlo << ' ' <<
/*7*/			ftildennlo << ' ' <<
/*8*/			deltaf1 << ' ' <<
/*9*/			deltaf2 << ' ' <<
/*10*/			deltaf3 << ' ' <<
/*11*/			deltafnlo << ' ' <<
/*12*/			deltafnnlo << '\n';
			
	}

	delete pdf;
}
