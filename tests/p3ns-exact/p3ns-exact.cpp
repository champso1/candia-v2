#include <filesystem>
using namespace std;
namespace fs = filesystem;

#include "Candia-v2/Candia.hpp"
using namespace Candia2;
using out_type = std::vector<ArrayGrid>;

static constexpr char const* DATAFILEDIR = "data";

static void usage();
static void outputData(std::vector<ArrayGrid> const& F, Grid const& grid, std::string filename);
static std::vector<ArrayGrid> calculate_ratios(
	std::vector<ArrayGrid> const& exact,
	std::vector<ArrayGrid> const& approx_central,
	std::vector<ArrayGrid> const& approx_imod1,
	std::vector<ArrayGrid> const& approx_imod2,
	int type);

int main()
{
    auto run = [&](std::vector<std::pair<ExprName, P3ApproxType>> const& approx_exprs) {
		const uint order = 3;
		const uint iterations = 12;
		const uint trunc_idx = 10;
		const double mur2_muf2 = 1.0;
		const double Qf = 100.0;

		auto& log_options = getLogOptions();
		log_options.verbosity = LOG_INFO;
	
		vector<double> xtab;
		Grid grid({1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0});

		LesHouchesDistribution dist(Qf);
		AlphaS alphas(order, dist.Q0(), Qf, dist.alpha0(), mur2_muf2);
		alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
		// alphas.setFFNS(4);

		DGLAPSolver solver(order, grid, alphas, Qf, iterations, trunc_idx, dist, mur2_muf2);
		auto& dglap_options = solver.getOptions();
		dglap_options.disable_heavy_flavor_matching = false;
		dglap_options.use_nnlo_matching_conditions_at_n3lo = false;
		dglap_options.use_n3lo_heavyquark_asymmetry = true;
		dglap_options.use_fortran_n3lo_splitfuncs = false;

		solver.setP3ApproximationTypes(approx_exprs);

		std::vector<ArrayGrid> F = solver.evolveTrunc();
		return F;
	};

	vector<double> xtab;
	Grid grid({1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0});

	std::vector<ArrayGrid> p3nsm_exact = run({std::make_pair(ExprName::P3nsm, P3ApproxType::Exact)});
	std::vector<ArrayGrid> p3nsm_approx_central = run({});
	std::vector<ArrayGrid> p3nsm_approx_imod1 = run({std::make_pair(ExprName::P3nsm, P3ApproxType::Imod1)});
	std::vector<ArrayGrid> p3nsm_approx_imod2 = run({std::make_pair(ExprName::P3nsm, P3ApproxType::Imod2)});
	std::vector<ArrayGrid> p3nsm_ratios_qns = calculate_ratios(p3nsm_exact, p3nsm_approx_central, p3nsm_approx_imod1, p3nsm_approx_imod2, 0);
	outputData(p3nsm_ratios_qns, grid, "p3nsm-exact-vs-approx-qns.dat");
	std::vector<ArrayGrid> p3nsm_ratios_lm = calculate_ratios(p3nsm_exact, p3nsm_approx_central, p3nsm_approx_imod1, p3nsm_approx_imod2, 2);
	outputData(p3nsm_ratios_lm, grid, "p3nsm-exact-vs-approx-lm.dat");

	std::vector<ArrayGrid> p3nsp_exact = run({std::make_pair(ExprName::P3nsp, P3ApproxType::Exact)});
	std::vector<ArrayGrid> p3nsp_approx_central = run({});
	std::vector<ArrayGrid> p3nsp_approx_imod1 = run({std::make_pair(ExprName::P3nsp, P3ApproxType::Imod1)});
	std::vector<ArrayGrid> p3nsp_approx_imod2 = run({std::make_pair(ExprName::P3nsp, P3ApproxType::Imod2)});
    std::vector<ArrayGrid> p3nsp_ratios_qns = calculate_ratios(p3nsp_exact, p3nsp_approx_central, p3nsp_approx_imod1, p3nsp_approx_imod2, 1);
	outputData(p3nsp_ratios_qns, grid, "p3nsp-exact-vs-approx-qns.dat");
	std::vector<ArrayGrid> p3nsp_ratios_lm = calculate_ratios(p3nsp_exact, p3nsp_approx_central, p3nsp_approx_imod1, p3nsp_approx_imod2, 2);
	outputData(p3nsp_ratios_lm, grid, "p3nsp-exact-vs-approx-lm.dat");
}

static void outputData(
	std::vector<ArrayGrid> const& F, Grid const& grid, std::string filename)
{
	fs::path datafile_path(filename);
	ofstream outfile(datafile_path);

	// print them out
	for (uint k=0; k<grid.size(); k++){
		outfile << setw(15) << setprecision(8) << grid[k] << ' ';
		outfile << setprecision(std::numeric_limits<double>::max_digits10);	
		for (uint j=0; j<F.size(); ++j)
			outfile << F[j][k] << ' ';
		outfile << '\n';
	}
}

static std::vector<ArrayGrid> calculate_ratios(
	std::vector<ArrayGrid> const& exact,
	std::vector<ArrayGrid> const& approx_central,
	std::vector<ArrayGrid> const& approx_imod1,
	std::vector<ArrayGrid> const& approx_imod2,
	int type)
{
	auto calc_qnsplus1d = [](std::vector<ArrayGrid> const& A, uint k, uint j){
		return (A[1][k] + A[1+6][k]) - (A[j][k] + A[j+6][k]);
	};
	auto calc_qnsminus1d = [](std::vector<ArrayGrid> const& A, uint k, uint j){
		return (A[1][k] - A[1+6][k]) - (A[j][k] - A[j+6][k]);
	};	
	auto calc_lm = [](std::vector<ArrayGrid> const& A, uint k, uint j){
		return (A[2+6][k] - A[1+6][k]);
	};
	
	auto size = exact[0].size();
	std::vector<ArrayGrid> ratios(4, ArrayGrid(exact[0].size()));
	for (uint k=0; k<size; ++k) {
		ratios[0][k] = 1.0;

		if (type == 0) {
			ratios[1][k] = calc_qnsminus1d(approx_central, k, 2)/calc_qnsminus1d(exact, k, 2);
			ratios[2][k] = calc_qnsminus1d(approx_imod1, k, 2)/calc_qnsminus1d(exact, k, 2);
			ratios[3][k] = calc_qnsminus1d(approx_imod2, k, 2)/calc_qnsminus1d(exact, k, 2);
		} else if (type == 1) {
			ratios[1][k] = calc_qnsplus1d(approx_central, k, 2)/calc_qnsplus1d(exact, k, 2);
			ratios[2][k] = calc_qnsplus1d(approx_imod1, k, 2)/calc_qnsplus1d(exact, k, 2);
			ratios[3][k] = calc_qnsplus1d(approx_imod2, k, 2)/calc_qnsplus1d(exact, k, 2);
		} else if (type == 2) {
			ratios[1][k] = calc_lm(approx_central, k, 2)/calc_lm(exact, k, 2);
			ratios[2][k] = calc_lm(approx_imod1, k, 2)/calc_lm(exact, k, 2);
			ratios[3][k] = calc_lm(approx_imod2, k, 2)/calc_lm(exact, k, 2);
		}
	}
	
	return ratios;
}
