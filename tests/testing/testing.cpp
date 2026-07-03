#include "Candia-v2/Candia.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"
#include <iomanip>
#include <limits>
using namespace Candia2;

#include "LHAPDF/LHAPDF.h"

constexpr double Q  = 10.0;
constexpr double Q2 = Q*Q;

int main()
{
	LHAPDF::PDF* pdf = LHAPDF::mkPDF("CT18NLO", 0);
	
	double as = pdf->alphasQ(10.0);
	double as2 = as*as;
	double as3 = as2*as;
	double mb = pdf->quarkMass(5);
	double mb2 = mb*mb;
	double L = std::log(Q2/mb2);
		
	OpMatElemN3LO::update(-L, 5);

	auto zero_func = [](double,double){ return 0.0; };
		
	// auto& p1qg = getExpression("P1qg");
	auto a1qg_reg_func = [as](double lm, double nf, double x) {
		auto trunced = ome::AQg_reg.truncate(1);
		return trunced(as, lm, nf, x); };
	OpMatElemCustom a1hg(a1qg_reg_func, zero_func, zero_func);

	auto a2hq_reg_func = [as](double lm, double nf, double x) {
		auto trunced = ome::AQqPS_reg.truncate(2);
		return trunced(as, lm, nf, x); };
	OpMatElemCustom a2hq(a2hq_reg_func, zero_func, zero_func);
		
	auto a2hg_reg_func = [as](double lm, double nf, double x) {
		auto trunced = ome::AQg_reg.truncate(2);
		return trunced(as, lm, nf, x); };
	OpMatElemCustom a2hg(a2hg_reg_func, zero_func, zero_func);

	auto a3hq_reg_func = [as](double lm, double nf, double x) {
		auto trunced = ome::AQqPS_reg.truncate(3);
		return trunced(as, lm, nf, x); };
	OpMatElemCustom a3hq(a3hq_reg_func, zero_func, zero_func);
		
	auto a3hg_reg_func = [as](double lm, double nf, double x) {
		auto trunced = ome::AQg_reg.truncate(3);
		return trunced(as, lm, nf, x); };
	OpMatElemCustom a3hg(a3hg_reg_func, zero_func, zero_func);

	Grid grid({1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0});

	ArrayGrid g(grid.size()-1), b(grid.size()-1);
	for (uint k=0; k<grid.size()-1; ++k) {
		g[k] = pdf->xfxQ(21, grid[k], Q);
		b[k] = pdf->xfxQ(5, grid[k], Q);
	}

	ArrayGrid sigma(grid.size()-1);
	for (uint k=0; k<grid.size()-1; ++k) {
		sigma[k] = 0.0;
		for (uint i=1; i<=5; ++i) {
			sigma[k] += pdf->xfxQ(i, grid[k], Q) + pdf->xfxQ(-i, grid[k], Q);
		}
	}

	std::vector<ArrayGrid> subpdfs(10, ArrayGrid(grid.size()-1));
	for (uint k=0; k<grid.size()-1; ++k) {
		double ftilde1 = as*grid.convolution(g, a1hg, k);
		double ftilde2 = as2*(
			grid.convolution(sigma, a2hq, k) +
			grid.convolution(g, a2hg, k)
		);
		double ftilde3 = as3*(
			grid.convolution(sigma, a3hq, k) +
			grid.convolution(g, a3hg, k)
		);

		subpdfs[0][k] = ftilde1;
		subpdfs[1][k] = ftilde2;
		subpdfs[2][k] = ftilde3;
			
		subpdfs[3][k] = ftilde1 + ftilde2;
		subpdfs[4][k] = ftilde1 + ftilde2 + ftilde3;
			
		subpdfs[5][k] = std::abs(b[k] - ftilde1);
		subpdfs[6][k] = std::abs(b[k] - ftilde2);
		subpdfs[7][k] = std::abs(b[k] - ftilde3);
		subpdfs[8][k]  = std::abs(b[k] - subpdfs[3][k]);
		subpdfs[9][k]  = std::abs(b[k] - subpdfs[4][k]);
	}

	std::ofstream outfile("out.dat");
	outfile << std::scientific;
	for (uint k=0; k<grid.size()-1; ++k) {
		outfile << grid[k] << ' ';
		outfile << b[k] << ' ';
		for (auto const& subpdf : subpdfs)
			outfile << subpdf[k] << ' ';
		outfile << '\n';
	}
}
	
