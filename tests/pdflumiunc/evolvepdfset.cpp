#include <filesystem>
using namespace std;
namespace fs = filesystem;

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/LHAPDFDistribution.hpp"
#include "Candia-v2/LHAPDFGrid.hpp"
using namespace Candia2;

int main()
{
	double q0 = 1.3;
	double qf = 14.0e3;
	double log_q0 = std::log10(q0);
	double log_qf = std::log10(qf);
	uint num = 30;
	double dnum = static_cast<double>(num);
    auto qvals_view =
		std::views::iota(uint{0}, num)
		| std::views::transform([&](uint i){ return std::pow(10.0, log_q0+(log_qf-log_q0)*static_cast<double>(i)/(num-1)); });
	std::vector<double> qvals(qvals_view.begin(), qvals_view.end());
	
	getLogOptions().verbosity = LOG_DEBUG;
	LHAPDFDistribution dist(make_lhapdf_pdf("CT25NNLO"), q0, qf);
	Grid grid({1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0});
	DGLAPOptions dglap_options{};
	dglap_options.use_truncated_nonsinglet_sol = true;
	dglap_options.disable_heavy_flavor_matching = false;
	dglap_options.use_nnlo_matching_conditions_at_n3lo = false;
	dglap_options.use_n3lo_heavyquark_asymmetry = true;
	dglap_options.use_fortran_n3lo_splitfuncs = false;

	std::string pdfname("CT25aN3LO");
	LHAPDFGrid lhapdfgrid(
		pdfname, fs::current_path()/"infofile.in",
		dist, grid,
		3, 15, 10, 1);
	
	lhapdfgrid.evolve(qvals, dglap_options);

	getLogOptions().verbosity = LOG_INFO;
	log(LOG_INFO, "evolvepdfset.cpp", "finished running the however many evolutions");
	lhapdfgrid.write();

	fs::path pdfdir = fs::current_path()/pdfname;
	fs::path lhapdfdir("/home/champson/.local/share/LHAPDF");
	fs::path pdfdir_dest = lhapdfdir/pdfname;
	if (!fs::exists(pdfdir_dest)) {
		if (!fs::create_directory(pdfdir_dest))
			log(LOG_ERROR, "evolvepdfset.cpp", "failed to create lhapdf output dir (in ~/.local/share/LHAPDF");
	}
	
	log(LOG_INFO, "evolvepdfset.cpp", "{} --> {}", pdfdir.string(), pdfdir_dest.string());
    const auto copyoptions = fs::copy_options::recursive | fs::copy_options::overwrite_existing;
	fs::copy(pdfdir, pdfdir_dest, copyoptions);
}
