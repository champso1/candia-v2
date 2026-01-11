#include "Candia-v2/Distribution.hpp"

namespace Candia2
{
    LHAPDFDistribution::LHAPDFDistribution(lhapdf_pdf_ptr_type lhapdf_pdf, double Q0)
        : Distribution(), _pdf(std::move(lhapdf_pdf)), _pids{lhapdf_pdf->flavors()}
    {
        _Q0 = Q0;
        for (uint i=0; i<_masses.size(); ++i) {
            double mass = _pdf->quarkMass(int id)
        }
    }

    double LHAPDFDistribution::xuv(double x) const
	{
		return 5.1072*std::pow(x, 0.8)*std::pow(1.0-x, 3.0);
	}

	double LHAPDFDistribution::xdv(double x) const
	{
		return 3.06432*std::pow(x, 0.8)*std::pow(1.0-x, 4.0);
	}

	double LHAPDFDistribution::xg(double x) const
	{
		return 1.7*std::pow(x, -0.1)*std::pow(1.0-x, 5.0);
	}

	double LHAPDFDistribution::xdb(double x) const
	{
		return 0.1939875*std::pow(x, -0.1)*std::pow(1.0-x, 6.0);
	}

	double LHAPDFDistribution::xub(double x) const
	{
		return xdb(x)*(1.0-x);
	}

	double LHAPDFDistribution::xs(double x) const
	{
		return 0.2*(xub(x) + xdb(x));
	}
}