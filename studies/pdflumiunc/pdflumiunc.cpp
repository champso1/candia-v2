#include "Candia-v2/Common.hpp"
#include "Candia-v2/Grid.hpp"
using namespace Candia2;

#include "LHAPDF/LHAPDF.h"

#include <iostream>
#include <fstream>
#include <filesystem>
#include <memory>
namespace fs = std::filesystem;

static constexpr double ECM = 13.6e3*13.6e3; // s

using unc_type = std::tuple<double,double,double,double>;
template <typename T>
using setval_type = std::vector<std::vector<T>>;

static setval_type<double> getGluonFromCandiaPDF(
	LHAPDF::PDF* pdf,
	Grid const& g,
	std::vector<double> const& q2vals);
static setval_type<std::vector<double>> getGluonFromPDF(
	LHAPDF::PDFSet const& s,
	Grid const& g,
	std::vector<double> const& q2vals);
static unc_type computeErrorForCalcedVals(
	LHAPDF::PDFSet const& set,
	std::vector<double> const& vals);


int main()
{
	LHAPDF::setVerbosity(0);
	// CT25aN3LO is written into this directory by evolvepdfset
	LHAPDF::pathsPrepend(fs::current_path().string());
	getLogOptions().verbosity = LOG_INFO;
	
	Grid grid({1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9});

	double q0 = 1.3;
	double qf = 13.59e3;
	double log_q0 = std::log10(q0);
	double log_qf = std::log10(qf);
	uint num = 200;
	// Q^2 values from q0^2 to qf^2, logarithmically spaced
    auto q2vals_view =
		std::views::iota(uint{0}, num)
		| std::views::transform([&](uint i){ return std::pow(std::pow(10.0, log_q0+(log_qf-log_q0)*static_cast<double>(i)/(num-1)), 2); });
	std::vector<double> q2vals(q2vals_view.begin(), q2vals_view.end());

	std::unique_ptr<LHAPDF::PDF> ct25an3lo(LHAPDF::mkPDF("CT25aN3LO", 0));
	LHAPDF::PDFSet msht("MSHT20an3lo_as118");
	LHAPDF::PDFSet nnpdf("NNPDF40_an3lo_as_01180_mhou_hessian");

    setval_type<double> ct25an3lo_xg = getGluonFromCandiaPDF(ct25an3lo.get(), grid, q2vals);
	setval_type<std::vector<double>> msht_xg = getGluonFromPDF(msht, grid, q2vals);
	setval_type<std::vector<double>> nnpdf_xg = getGluonFromPDF(nnpdf, grid, q2vals);

	std::vector<double> ct25an3lo_gg_lumis(q2vals.size());
	setval_type<double> msht_gg_lumis(q2vals.size(), setval_type<double>::value_type(msht.size()));
	setval_type<double> nnpdf_gg_lumis(q2vals.size(), setval_type<double>::value_type(nnpdf.size()));
	
	for (uint qi=0; qi<q2vals.size(); ++qi) {
		double q2 = q2vals[qi];
		double tau = q2/ECM;
		double fac = 1.0/q2;

		std::span<double> ct25an3lo_vals(ct25an3lo_xg[qi]);
		ct25an3lo_gg_lumis[qi] = fac*grid.convolution(ct25an3lo_vals, ct25an3lo_vals, tau);

		std::vector<std::vector<double>>& msht_vals = msht_xg[qi];
		for (uint i=0; i<msht.size(); ++i) {
			std::span<double> msht_set_vals(msht_vals[i]);
			msht_gg_lumis[qi][i] = fac*grid.convolution(msht_set_vals, msht_set_vals, tau);
		}

		std::vector<std::vector<double>>& nnpdf_vals = nnpdf_xg[qi];
		for (uint i=0; i<nnpdf.size(); ++i) {
			std::span<double> nnpdf_set_vals(nnpdf_vals[i]);
			nnpdf_gg_lumis[qi][i] = fac*grid.convolution(nnpdf_set_vals, nnpdf_set_vals, tau);
		}
	}
	
	fs::path outfile_path("lumis.dat");
	std::ofstream outfile(outfile_path);
	for (uint qi=0; qi<q2vals.size(); ++qi) {
		double xaxis_val = std::sqrt(q2vals[qi]/ECM);

		unc_type msht_lumis = computeErrorForCalcedVals(msht, msht_gg_lumis[qi]);
		unc_type nnpdf_lumis = computeErrorForCalcedVals(nnpdf, nnpdf_gg_lumis[qi]);
		
		double baseline = std::get<0>(nnpdf_lumis);
		
		outfile << xaxis_val << ' ' << ct25an3lo_gg_lumis[qi]/baseline << ' ';
	    outfile
			<< std::get<0>(msht_lumis)/baseline << ' '
			<< std::get<1>(msht_lumis)/baseline << ' '
			<< std::get<2>(msht_lumis)/baseline << ' ';
		
	    outfile
			<< std::get<0>(nnpdf_lumis)/baseline << ' '
			<< std::get<1>(nnpdf_lumis)/baseline << ' '
			<< std::get<2>(nnpdf_lumis)/baseline << ' ';

		outfile << std::endl;
	}
	
}


static
setval_type<double>
getGluonFromCandiaPDF(
	LHAPDF::PDF* pdf,
	Grid const& g,
	std::vector<double> const& q2vals)
{
    setval_type<double> xg(q2vals.size(), setval_type<double>::value_type(g.size(), 0.0));
    for (uint qi=0; qi<q2vals.size(); ++qi) {
		for (uint xi=0; xi<g.size(); ++xi)
			xg[qi][xi] = pdf->xfxQ2(21,g[xi],q2vals[qi]);
	}
	return xg;
}

static
setval_type<std::vector<double>>
getGluonFromPDF(
	LHAPDF::PDFSet const& set,
	Grid const& g,
	std::vector<double> const& q2vals)
{
	const std::vector<LHAPDF::PDF*> pdfs = set.mkPDFs();
    setval_type<std::vector<double>> xg(
		q2vals.size(), std::vector<std::vector<double>>(
			set.size(), std::vector<double>(g.size(), 0.0)));
	for (uint qi=0; qi<q2vals.size(); ++qi) {
		for (uint i=0; i<set.size(); ++i) {
			for (uint xi=0; xi<g.size(); ++xi)
				xg[qi][i][xi] = pdfs[i]->xfxQ2(21,g[xi],q2vals[qi]);
		}
	}
	for (LHAPDF::PDF* pdf : pdfs)
		delete pdf;
	return xg;
}


static unc_type computeErrorForCalcedVals(
	LHAPDF::PDFSet const& set,
	std::vector<double> const& vals)
{
	const LHAPDF::PDFUncertainty err = set.uncertainty(vals, -1.0);
	return unc_type{err.central, err.errplus, err.errminus, err.errsymm};
}
