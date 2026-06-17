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

	LesHouchesDistribution dist(100.0);
	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	GridFillerLogLinQuad grid_filler(1e-5, 100, 50, 25);
	Grid grid(xtab, grid_filler, {});
	DGLAPOptions dglap_options{};
	dglap_options.use_truncated_nonsinglet_sol = false;
	dglap_options.disable_heavy_flavor_matching = false;
	dglap_options.use_nnlo_matching_conditions_at_n3lo = false;
	dglap_options.use_n3lo_heavyquark_asymmetry = true;
	dglap_options.use_fortran_n3lo_splitfuncs = false;
	
	LHAPDFGrid lhapdfgrid(
		"testpdf", fs::current_path()/"infofile.in",
		dist, grid,
		3, 10, 10, 1);
	
	lhapdfgrid.evolve(std::numbers::sqrt2, 110.0, 5.0, dglap_options);

	getLogOptions().verbosity = LOG_INFO;
	lhapdfgrid.write();

	fs::path testpdfdir = fs::current_path()/"testpdf";
	fs::path testpdfdir_dest = fs::path("/home/champson/.local/share/LHAPDF/testpdf");
	log(LOG_INFO, "lhapdf-datafile.cpp", "{} --> {}", testpdfdir.string(), testpdfdir_dest.string());
    const auto copyoptions = fs::copy_options::recursive | fs::copy_options::overwrite_existing;
	fs::copy(testpdfdir, testpdfdir_dest, copyoptions);
}
