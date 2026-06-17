#include "Candia-v2/LHAPDFGrid.hpp"

#include "LHAPDF/LHAPDF.h"

#include <vector>
#include <fstream>
#include <format>
#include <ranges>

int main()
{
	LHAPDF::PDF* pdf = LHAPDF::mkPDF("testpdf", 0);
	
	double xmin = 1e-5;
	double xmax = 0.99;
	double logxmin = std::log10(xmin);
	double logxmax = std::log10(xmax);
	uint nx = 200;
	double delta = (logxmax-logxmin)/nx;
	auto xvals =
		std::views::iota(uint{0}, nx)
		| std::views::transform(
			[&](uint _i) {
				double i = static_cast<double>(_i);
				return std::pow(10.0, logxmin + delta*i);
			});
	
	std::vector<double> qs{10.0, 100.0};
	for (double q : qs) {
		std::ofstream outfile(std::format("out-q{}.dat", q));
		for (double x : xvals) {
			outfile <<
				x << ' ' <<
				pdf->xfxQ(5, x, q) << ' ' <<
				pdf->xfxQ(Candia2::LHAPDFGrid::FTILDE1, x, q) << ' ' <<
				pdf->xfxQ(Candia2::LHAPDFGrid::FTILDE2, x, q) << ' ' <<
				pdf->xfxQ(Candia2::LHAPDFGrid::FTILDE3, x, q) << ' ' <<
				pdf->xfxQ(Candia2::LHAPDFGrid::FTILDENNLO, x, q) << ' ' <<
				pdf->xfxQ(Candia2::LHAPDFGrid::FTILDEN3LO, x, q) << ' ' <<
				pdf->xfxQ(Candia2::LHAPDFGrid::DELTAF1, x, q) << ' ' <<
				pdf->xfxQ(Candia2::LHAPDFGrid::DELTAF2, x, q) << ' ' <<
				pdf->xfxQ(Candia2::LHAPDFGrid::DELTAFNNLO, x, q) << ' ' <<
				pdf->xfxQ(Candia2::LHAPDFGrid::DELTAFN3LO, x, q) << '\n';
		}
	}

	delete pdf;
}
