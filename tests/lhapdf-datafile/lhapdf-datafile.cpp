#include "Candia-v2/Candia.hpp"
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
		1, 10, 10, 1.0);
	GridOptions grid_options{};
	grid_options.use_alt_mapping = true;
	grid_options.use_gsl_conv_routine = false;
	grid_options.use_gsl_interp_routine = true;
	
	DGLAPOptions dglap_options{};
	dglap_options.use_truncated_nonsinglet_sol = false;
	dglap_options.disable_heavy_flavor_matching = false;
	dglap_options.use_nnlo_matching_conditions_at_n3lo = false;
	dglap_options.use_n3lo_heavyquark_asymmetry = true;
	dglap_options.use_fortran_n3lo_splitfuncs = false;
	dglap_options.cache_exprs = false;
	
	solver.evolve(std::numbers::sqrt2, 100, 5.0, grid_options, dglap_options);
	solver.write();
}
