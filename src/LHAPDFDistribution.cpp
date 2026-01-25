#include "Candia-v2/Distribution.hpp"

namespace Candia2
{
    LHAPDFDistribution::LHAPDFDistribution(lhapdf_pdf_ptr_type lhapdf_pdf, double Q0)
        : Distribution(),
		_pdf(std::move(lhapdf_pdf)), _pids{lhapdf_pdf->flavors()}
    {
        _Q0 = Q0;
        for (uint i=0; i<6; ++i) {
            double mass = _pdf->quarkMass(i+1);
			if (mass >= Q0)
				_masses[i] = mass;
        }
    }

    double LHAPDFDistribution::xuv(double x) const
	{
		return _pdf->xfxQ(1, x, _Q0) - _pdf->xfxQ(-1, x, _Q0);
	}

	double LHAPDFDistribution::xdv(double x) const
	{
		return _pdf->xfxQ(2, x, _Q0) - _pdf->xfxQ(-2, x, _Q0);
	}

	double LHAPDFDistribution::xg(double x) const
	{
		return _pdf->xfxQ(21, x, _Q0);
	}

	double LHAPDFDistribution::xdb(double x) const
	{
		return _pdf->xfxQ(-2, x, _Q0);
	}

	double LHAPDFDistribution::xub(double x) const
	{
		return _pdf->xfxQ(-1, x, _Q0);
	}

	double LHAPDFDistribution::xs(double x) const
	{
		return _pdf->xfxQ(3, x, _Q0);
	}
}