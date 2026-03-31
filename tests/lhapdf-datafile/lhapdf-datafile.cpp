#include "Candia-v2/Candia.hpp"
using namespace Candia2;

#include <filesystem>
using namespace std;

int main()
{
	const uint order = 3;
	const uint iterations = 10;
	const uint trunc_idx = 10;
	const double Qf = 100.0;
	const double mur2_muf2 = 1.0;

	auto logopts = getLogOptions();
	logopts.show_debug_messages = true;
	setLogOptions(logopts);
	
	vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab,
		make_grid_filler<GridFillerLogLinQuad>(1e-5, 101, 51, 26),
		{ .default_gauss_points=70, .split_interval = true});
	auto& grid_options = grid.getOptions();
	grid_options.use_alt_mapping = true;
	grid_options.use_gsl_conv_routine = false;
	grid_options.use_gsl_interp_routine = true;

	LesHouchesDistribution dist{};
	AlphaS alphas(order, dist.Q0(), Qf, dist.alpha0(), mur2_muf2);
	alphas.setVFNS(dist.masses(), dist.nfi());

	DGLAPSolver solver(order, grid, alphas, Qf, iterations, trunc_idx, dist, mur2_muf2);
	auto& dglap_options = solver.getOptions();
	dglap_options.use_truncated_nonsinglet_sol = true;
	dglap_options.disable_heavy_flavor_matching = false;
	dglap_options.use_nnlo_matching_conditions_at_n3lo = false;
	dglap_options.use_n3lo_heavyquark_asymmetry = true;
	dglap_options.use_fortran_nnlo_splitfuncs = false;
	dglap_options.use_fortran_n3lo_splitfuncs = true;
	dglap_options.cache_exprs = true;

	std::filesystem::path infofile_in_path("infofile.in");
	solver.generateLHAPDFGrid("testpdf", infofile_in_path);
}
