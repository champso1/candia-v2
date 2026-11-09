// Bottom-quark subtraction PDFs with CT18 NLO and NNLO inputs,
// as in the first three pages of the working note by Keping
// We first fix the scale of the PDFs to 10 and 100,
// then print the b quark, LO and NLO subtraction PDFs, and their residuals
// for CT18NLO and NNLO -> 4 plots (page 1)
// then for the two PDF sets, we fix x at 1e-4, 1e-2, and 0.1,
// and perform a logarithmic scan up to 10TeV -> 6 plots (page 2)
// page 1 is done again but with just a single fixed energy, 5Gev,
// right above the bottom quark mass -> 2 plots (page 3)
// we also spit out two plots containing also the NNLO subtraction PDFs
// and their residuals as 10 and 100 GeV for the CPC paper



#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/LHAPDFDistribution.hpp"
#include "Candia-v2/LHAPDFGrid.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <format>
#include <iostream>
#include <string>
#include <vector>
#include <ranges>

using namespace Candia2;

namespace
{
	// x f_b, x ftilde^LO, x ftilde^NLO, x delta f^LO, x delta f^NLO
	constexpr uint NUM_QUANTITIES = 5;

	const std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};

	struct PDFSet
	{
		std::string name;
		uint order;
	};

	// values[scale][quantity][grid point]
	using MemberValues = std::vector<std::array<std::vector<double>, NUM_QUANTITIES>>;

	std::vector<double> scales(double mb)
	{
		constexpr uint num_scan = 30;
		constexpr double q_max = 1.0e4;
		std::vector<double> qs{mb, 5.0, 10.0, 100.0};
		for (uint i=1; i<num_scan; ++i)
			qs.push_back(mb*std::pow(q_max/mb, static_cast<double>(i)/(num_scan-1)));
		std::ranges::sort(qs);
		return qs;
	}

	void evolveMember(PDFSet const& set, LHAPDFDistribution& dist, Grid grid,
		std::vector<double> const& qs, std::vector<std::array<std::vector<double>, NUM_QUANTITIES>>& out)
	{
		constexpr uint iterations = 15;
		constexpr uint trunc_idx = 10;
		const double q0 = dist.Q0(); // the set's qMin(), see main()

		out.assign(qs.size(), {});
		for (uint iq=0; iq<qs.size(); ++iq) {
			const double q = qs[iq];
			dist.setup(q0, q);
			AlphaS alphas(set.order, q0, q, dist.alpha0(), 1.0);
			// CT18 has five flavours; at m_b itself this gives the values just above the threshold
			alphas.setVFNS(dist.masses(), dist.nfi(), 5);

			DGLAPSolver solver(set.order, grid, alphas, q, iterations, trunc_idx, dist, 1.0);
			auto const& F = solver.evolve();
			auto const sub = solver.calculateSubtractionPDFs();

			for (auto& values : out[iq])
				values.resize(grid.size());
			for (uint k=0; k<grid.size(); ++k) {
				const double fb = 0.5*(F[DIST_B][k] + F[DIST_B+6][k]);
				const double ft_lo = sub[LHAPDFGrid::B_FTILDE1][k];
				const double ft_nlo = sub[LHAPDFGrid::B_FTILDENLO][k];
				out[iq][0][k] = fb;
				out[iq][1][k] = ft_lo;
				out[iq][2][k] = ft_nlo;
				out[iq][3][k] = fb - ft_lo;
				out[iq][4][k] = fb - ft_nlo;
			}
		}
	}

	void writeSubtractionPDFs(double Q, uint heavy, std::filesystem::path const& filepath)
	{
		constexpr uint order = 3;
		constexpr uint iterations = 20;
		constexpr uint trunc_idx = 15;
		constexpr double mur2_muf2 = 1;
		
		Grid grid(xtab);
		LesHouchesDistribution dist(Q);
		AlphaS alphas(order, dist.Q0(), Q, dist.alpha0(), mur2_muf2);
		alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
		DGLAPSolver solver(order, grid, alphas, Q, iterations, trunc_idx, dist, mur2_muf2);
		auto const& F = solver.evolve();
		std::vector<ArrayGrid> subtraction_pdfs = solver.calculateSubtractionPDFs();
		const uint first = (heavy == DIST_C) ? 0 : 10; // charm entries are 0-9, bottom 10-19

		std::ofstream outfile(filepath);
		outfile << "# subtraction PDFs at Q = " << Q << " GeV: x f f~1 f~2 f~3 f~NLO f~NNLO f-f~1 f-f~2 f-f~3 f-f~NLO f-f~N3LO\n";
		outfile << std::scientific << std::setprecision(6);
		for (uint k=0; k<grid.size(); ++k) {
			outfile << grid[k] << ' ' << F[heavy][k];
			for (uint i=first; i<first+10; ++i)
				outfile << ' ' << subtraction_pdfs[i][k];
			outfile << '\n';
		}
		if (!outfile)
			log(LOG_ERROR, "subtraction-pdfs.cpp", "failed to write '{}'", filepath.string());
	}
}

int main()
{
	const std::vector<PDFSet> sets{{"CT18NLO", 1}, {"CT18NNLO", 2}};

	getLogOptions().verbosity = LOG_WARNING;
	LHAPDF::setVerbosity(0);
	std::filesystem::create_directories("data");

	Grid grid(xtab);

	for (uint i=0; i<sets.size(); ++i) {
		auto const& set = sets[i];
		auto pdf = make_lhapdf_pdf(set.name);
		const double q0 = pdf->qMin();
		std::unique_ptr<LHAPDFDistribution> dist = std::make_unique<LHAPDFDistribution>(std::move(pdf), q0, q0);
		const double mb = dist->masses(DIST_B);
		const std::vector<double> qs = scales(mb);

		std::cout << std::format("{}: evolving to {} scales at order {}\n", set.name, qs.size(), set.order);
		MemberValues values(qs.size());
		evolveMember(set, *dist, grid, qs, values);

		std::ofstream out("data/" + set.name + ".dat");
		out << std::format("# {} evolved by Candia-v2 at order {} from Q0 = {} GeV\n",
			set.name, set.order, dist->Q0());
		out << std::format("# m_b = {} GeV; Q = m_b is just above the threshold; xb_lhapdf is the central member's own grid value\n",
			mb);
		out << std::format("# {:>12} {:>14} {:>14} {:>14} {:>14} {:>14} {:>14} {:>14}\n",
			"Q", "x", "xfb", "xft_LO", "xft_NLO", "xdf_LO", "xdf_NLO", "xb_lhapdf");
		auto const& points = grid.points();
		for (uint iq=0; iq<qs.size(); ++iq) {
			for (uint k=0; k<grid.size(); ++k) {
				out << std::format("{:14.8e} {:14.8e}", qs[iq], points[k]);
				for (uint i=0; i<NUM_QUANTITIES; ++i)
					out << std::format(" {:14.7e} ", values[iq][i][k]);
				
				out << std::format(" {:14.7e}\n", dist->pdf().xfxQ(5, points[k], qs[iq]));
			}
		}
	}

	getLogOptions().verbosity = LOG_INFO;
	log(LOG_INFO, "subtraction-pdfs.cpp", "finished pdfset generation, writing specific grids for CPC plots");
	getLogOptions().verbosity = LOG_WARNING;
	
	writeSubtractionPDFs(4, DIST_C, "data/out-c-4gev.dat");
	writeSubtractionPDFs(10, DIST_B, "data/out-b-10gev.dat");
}
