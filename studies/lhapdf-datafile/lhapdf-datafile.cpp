// lhapdf-datafile.cpp
//
// Writes an N3LO LHAPDF grid of the Les Houches toy PDFs, which includes the subtraction PDFs
// (PDG ids 9000+), and tabulates the charm subtraction PDFs at Q = 4 GeV and the bottom ones at
// Q = 10 GeV for plot.gplt.

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/LHAPDFGrid.hpp"
using namespace Candia2;

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace
{
	constexpr uint order = 3;
	constexpr uint iterations = 10;
	constexpr uint trunc_idx = 10;
	constexpr double mur2_muf2 = 1.0;
	const std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};

	/**
	 *  @brief writes the subtraction PDFs of the heavy quark @a heavy (DIST_C or DIST_B) at the energy @a Q
	 *  columns: x, f, f~(1), f~(2), f~(3), f~(1)+f~(2), f~(1)+f~(2)+f~(3), and f minus each of the five f~
	 */
	void writeSubtractionPDFs(double Q, uint heavy, std::string const& filename)
	{
		Grid grid(xtab);
		LesHouchesDistribution dist(Q);
		AlphaS alphas(order, dist.Q0(), Q, dist.alpha0(), mur2_muf2);
		alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
		DGLAPSolver solver(order, grid, alphas, Q, iterations, trunc_idx, dist, mur2_muf2);
		auto const& F = solver.evolveTrunc();
		std::vector<ArrayGrid> subtraction_pdfs = solver.calculateSubtractionPDFs();
		const uint first = (heavy == DIST_C) ? 0 : 10; // charm entries are 0-9, bottom 10-19

		std::ofstream outfile(filename);
		outfile << "# subtraction PDFs at Q = " << Q << " GeV: x f f~1 f~2 f~3 f~12 f~123 f-f~1 f-f~2 f-f~3 f-f~12 f-f~123\n";
		outfile << std::scientific << std::setprecision(6);
		for (uint k=0; k<grid.size(); ++k) {
			outfile << grid[k] << ' ' << F[heavy][k];
			for (uint i=first; i<first+10; ++i)
				outfile << ' ' << subtraction_pdfs[i][k];
			outfile << '\n';
		}
		if (!outfile)
			log(LOG_ERROR, "lhapdf-datafile.cpp", "failed to write '{}'", filename);
	}
}

int main()
{
	getLogOptions().verbosity = LOG_WARNING;

	// the LHAPDF grid, with Q knots logarithmically spaced from Q0 to 100 GeV;
	// LHAPDFGrid splits it at m_b itself
	constexpr double Qmax = 100.0;
	{
		LesHouchesDistribution dist(Qmax);
		Grid grid(xtab);
		constexpr uint num_knots = 20;
		std::vector<double> qvals{};
		for (uint i=0; i<num_knots; ++i)
			qvals.push_back(dist.Q0()*std::pow(Qmax/dist.Q0(), static_cast<double>(i)/(num_knots-1)));

		LHAPDFGrid lhapdfgrid("testpdf", "infofile.in", dist, grid, order, iterations, trunc_idx, mur2_muf2);
		lhapdfgrid.evolveTrunc(qvals, {});
		lhapdfgrid.write();
	}
	std::cout << "wrote the LHAPDF set ./testpdf; add this directory to LHAPDF_DATA_PATH to use it\n";

	writeSubtractionPDFs(4.0, DIST_C, "out-c.dat");
	writeSubtractionPDFs(10.0, DIST_B, "out-b.dat");
	std::cout << "wrote out-c.dat and out-b.dat; plot them with 'gnuplot plot.gplt'\n";
}
