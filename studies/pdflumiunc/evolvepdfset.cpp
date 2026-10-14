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
    auto qvals_view =
		std::views::iota(uint{0}, num)
		| std::views::transform([&](uint i){ return std::pow(10.0, log_q0+(log_qf-log_q0)*static_cast<double>(i)/(num-1)); });
	std::vector<double> qvals(qvals_view.begin(), qvals_view.end());
	
	getLogOptions().verbosity = LOG_INFO;
	LHAPDFDistribution dist(make_lhapdf_pdf("CT25NNLO"), q0, qf);
	Grid grid({1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9});

	uint order = 3;
	uint iterations = 15;
	uint trunc_idx = 10;
	double mur2_muf2 = 1.0;

	std::string pdfname("CT25aN3LO");
	LHAPDFGrid lhapdfgrid(
		pdfname, fs::current_path()/"infofile.in",
		dist, grid,
		order, iterations, trunc_idx, mur2_muf2);
	
	lhapdfgrid.evolveTrunc(qvals, {});

	lhapdfgrid.write();

	// pdflumiunc, run from this directory, finds the set here; to use it elsewhere,
	// add this directory to LHAPDF_DATA_PATH or copy the set into an LHAPDF data directory
	log(LOG_INFO, "evolvepdfset.cpp", "wrote the LHAPDF set {}", (fs::current_path()/pdfname).string());
}
