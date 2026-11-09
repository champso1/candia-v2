// lhapdf-datafile.cpp
//
// Writes an N3LO LHAPDF grid of the Les Houches toy PDFs, which includes the subtraction PDFs
// (PDG ids 9000+), and tabulates the charm subtraction PDFs at Q = 4 GeV and the bottom ones at
// Q = 10 GeV for plot.gplt.

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/LHAPDFDistribution.hpp"
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
	constexpr uint order = 1;
	constexpr uint iterations = 20;
	constexpr uint trunc_idx = 15;
	constexpr double mur2_muf2 = 1.0;
	const std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};

	constexpr uint num_subpdfs_per_hq = 10; // 12: ftilde1-3, ftildeNLO,NNLO (5), and their residuals (5*2=10)

	/**
	 *  @brief writes the subtraction PDFs of the heavy quark @a heavy (DIST_C or DIST_B) at the energy @a Q
	 *  columns: x, f, f~(1), f~(2), f~(3), f~(LO), f~(NLO), f~(NNLO), and f minus each these 
	 */
	void writeSubtractionPDFs(double Q, uint heavy, std::filesystem::path const& filepath)
	{
		Grid grid(xtab);
		LesHouchesDistribution dist(Q);
		AlphaS alphas(order, dist.Q0(), Q, dist.alpha0(), mur2_muf2);
		alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
		DGLAPSolver solver(order, grid, alphas, Q, iterations, trunc_idx, dist, mur2_muf2);
		auto const& F = solver.evolve();
		std::vector<ArrayGrid> subtraction_pdfs = solver.calculateSubtractionPDFs();
		const uint first = (heavy == DIST_C) ? 0 : num_subpdfs_per_hq; // charm entries are 0-9, bottom 10-19

		std::ofstream outfile(filepath);
		outfile << "# subtraction PDFs at Q = " << Q << " GeV: x f f~1 f~2 f~3 f~NLO f~NNLO f-f~1 f-f~2 f-f~3 f-f~NLO f-f~N3LO\n";
		outfile << std::scientific << std::setprecision(6);
		for (uint k=0; k<grid.size(); ++k) {
			outfile << grid[k] << ' ' << F[heavy][k];
			for (uint i=first; i<first+num_subpdfs_per_hq; ++i)
				outfile << ' ' << subtraction_pdfs[i][k];
			outfile << '\n';
		}
		if (!outfile)
			log(LOG_ERROR, "lhapdf-datafile.cpp", "failed to write '{}'", filepath.string());
	}
}

int main()
{
	getLogOptions().verbosity = LOG_WARNING;

	// the LHAPDF grid, with Q knots logarithmically spaced from Q0 to 100 GeV;
	// LHAPDFGrid splits it at m_b itself
	constexpr double Qmax = 100.0;
	{
		// LesHouchesDistribution dist(Qmax);
		LHAPDFDistribution dist(make_lhapdf_pdf("CT18NNLO"), 1.295, Qmax);
		Grid grid(xtab);
		constexpr uint num_knots = 20;
		std::vector<double> qvals{};
		for (uint i=0; i<num_knots; ++i)
			qvals.push_back(dist.Q0()*std::pow(Qmax/dist.Q0(), static_cast<double>(i)/(num_knots-1)));

		LHAPDFGrid lhapdfgrid("testpdf-nlo", "infofile.in", dist, grid, order, iterations, trunc_idx, mur2_muf2);
		lhapdfgrid.evolve(qvals, {});
		lhapdfgrid.write();
	}

	
	writeSubtractionPDFs(2.0, DIST_C, "./data-nlo/out-c-2gev.dat");
	writeSubtractionPDFs(4.0, DIST_C, "./data-nlo/out-c-4gev.dat");
	writeSubtractionPDFs(10.0, DIST_B, "./data-nlo/out-b-10gev.dat");
}
