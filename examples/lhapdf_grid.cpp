#include "Candia-v2/Candia.hpp"
#include "Candia-v2/LHAPDFGrid.hpp"
#include "Candia-v2/LHAPDFDistribution.hpp"
using namespace Candia2;

#include <cmath>
#include <iostream>
#include <vector>

#include "LHAPDF/LHAPDF.h"

int main()
{
	constexpr double Qmax = 100.0;

	// create the standard LHAPDF toy model initial conditions
	// every evolution starts at the initial energy of the distribution, here Q0 = sqrt(2) GeV
	LesHouchesDistribution dist(Qmax);

	// the Q knots of the grid: logarithmically spaced from Q0 to Qmax.
	// LHAPDFGrid adds the heavy-quark thresholds (here m_b = 4.5 GeV) itself,
	// and splits the grid there so the PDFs can be discontinuous across them
	constexpr uint num_knots = 25;
	std::vector<double> qvals{};
	for (uint i=0; i<num_knots; ++i)
		qvals.push_back(dist.Q0()*std::pow(Qmax/dist.Q0(), static_cast<double>(i)/(num_knots-1)));

	// create the x-grid using the traditional xtab array
	Grid grid({
			1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9});

	/* evolution setup */
	/* see evolve_dglap.cpp for further description of these variables */
	uint order = 1;
	uint iterations = 15;
	uint trunc_idx = 15;
	double mur2_muf2 = 1.0;

	// create the LHAPDFGrid object called "testpdf"
	// 'infofile.in' only provides metadata (description, authors, ...);
	// everything that depends on the grid is generated
	LHAPDFGrid lhapdfgrid(
		"testpdf", "infofile.in",
		dist, grid,
		order, iterations, trunc_idx, mur2_muf2);
	// lhapdfgrid.evolve(qvals, {});   // exact
	lhapdfgrid.evolveTrunc(qvals, {}); // truncated
	lhapdfgrid.write(); // writes the PDF to disk

	// load the PDF back in
	// we can add the current directory to the search paths with this function
	// also remove the LHAPDF splash text
	LHAPDF::setVerbosity(0);
	LHAPDF::pathsPrepend(".");
	lhapdf_pdf_ptr_type testpdf = make_lhapdf_pdf("testpdf", 0);

	// print back the gluon at x=0.1, Q=75 GeV to see if it works
	std::cout << "testpdf value xg(x=0.1, Q=75) = " << testpdf->xfxQ(21, 0.1, 75.0) << std::endl;
	std::cout << "testpdf value alpha_s(Q=75)   = " << testpdf->alphasQ(75.0) << std::endl;
}
