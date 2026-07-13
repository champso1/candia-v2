#include "Candia-v2/Candia.hpp"
#include "Candia-v2/LHAPDFGrid.hpp"
using namespace Candia2;

#include <iostream>
#include <vector>

#include "LHAPDF/LHAPDF.h"

int main()
{
	// we want to evolve to 10GeV, 50GeV, and 100GeV
	std::vector<double> qvals{10.0, 50.0, 100.0};
	// create the standard LHAPDF toy model initial conditions
	// with the final energy being the last energy in the above list
	LesHouchesDistribution dist(qvals.back());
	// create the x-grid using the traditional xtab array
	Grid grid({
			1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9});

	/* evolution setup */
	/* see evolve_dglap.cpp for further description of these variables */
	uint order = 1;
	uint iterations = 10;
	uint trunc_idx = 10;
	double mur2_muf2 = 1.0;

	// create the LHAPDFGrid object called "testpdf"
	// using 'infofile.in' as the template info file for the LHAPDF PDF
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
	LHAPDF::paths().emplace_back(".");
	LHAPDF::PDF* testpdf = LHAPDF::mkPDF("testpdf", 0);

	// just print back the gluon at x=0.1,Q=75.0 to see if it works
	double x = 0.1;
	double Q2 = 75.0*75.0;
	std::cout << "testpdf value xg(x=0.1,Q2=75^2) = " << testpdf->xfxQ2(21, x, Q2) << std::endl;

	delete testpdf;
}
