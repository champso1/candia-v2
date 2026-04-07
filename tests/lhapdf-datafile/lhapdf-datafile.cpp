#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Grid.hpp"
using namespace Candia2;
using namespace std;

int main()
{
	auto& log_options = getLogOptions();
	log_options.show_debug_messages = true;
	log_options.show_thread_output = true;

	std::unique_ptr<Distribution> dist = std::make_unique<LesHouchesDistribution>();
	std::filesystem::path infofile_in_path("infofile.in");
	DGLAPSolverLHAPDF solver(
		"testpdf", infofile_in_path,
		std::move(dist),
		2, 6, 7, 1.0);
	solver.addSubtractionPDFs();
	
	DGLAPOptions dglap_options{};
	dglap_options.use_truncated_nonsinglet_sol = true;
	dglap_options.disable_heavy_flavor_matching = false;
	dglap_options.use_nnlo_matching_conditions_at_n3lo = false;
	dglap_options.use_n3lo_heavyquark_asymmetry = true;
	dglap_options.use_fortran_n3lo_splitfuncs = false;
	dglap_options.cache_exprs = false;

	vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
    auto grid_filler = make_grid_filler<GridFillerLogLinQuad>(1e-5, 101, 51, 26);
	GaussLegendreArgs gauleg_args{ .default_gauss_points=70, .split_interval = true};
	GridOptions grid_options{};
	grid_options.use_alt_mapping = true;
	grid_options.use_gsl_conv_routine = false;
	grid_options.use_gsl_interp_routine = true;
	
	solver.setGridInfo(xtab, std::move(grid_filler), gauleg_args, grid_options);
	
	solver.evolve(std::numbers::sqrt2, 100, 20.0, dglap_options);
	solver.write();
}
