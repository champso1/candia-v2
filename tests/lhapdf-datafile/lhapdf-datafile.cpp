#include "Candia-v2/Common.hpp"
#include "Candia-v2/LHAPDFGrid.hpp"
#include <LHAPDF/Config.h>
#include <LHAPDF/Factories.h>
#include <LHAPDF/LHAGlue.h>
#include <filesystem>
using namespace Candia2;

#include "LHAPDF/LHAPDF.h"
using namespace LHAPDF;

namespace fs = std::filesystem;

int main()
{
	getLogOptions().verbosity = LOG_WARNING;

	double qf = 100.0;
	LesHouchesDistribution dist(qf);
	std::vector<double> qvals{10.0, 50.0, 100.0};
	
	Grid grid({1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0});

	uint order = 3;
	uint iterations = 10;
	uint trunc_idx = 10;
	double mur2_muf2 = 1.0;
	
	LHAPDFGrid lhapdfgrid(
		"testpdf", "infofile.in",
		dist, grid,
		order, iterations, trunc_idx, mur2_muf2);

	// lhapdfgrid.evolve(qvals, {});
	lhapdfgrid.evolveTrunc(qvals, {});

	getLogOptions().verbosity = LOG_INFO;
	lhapdfgrid.write();

	fs::path testpdfdir = fs::current_path()/"testpdf";
	fs::path testpdfdir_dest = fs::path("/home/champson/.local/share/LHAPDF/testpdf");
	log(LOG_INFO, "lhapdf-datafile.cpp", "{} --> {}", testpdfdir.string(), testpdfdir_dest.string());
    const auto copyoptions = fs::copy_options::recursive | fs::copy_options::overwrite_existing;
	fs::copy(testpdfdir, testpdfdir_dest, copyoptions);
}
