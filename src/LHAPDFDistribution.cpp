#include "Candia-v2/LHAPDFDistribution.hpp"

namespace Candia2
{
	LHAPDFDistribution::LHAPDFDistribution(lhapdf_pdf_ptr_type lhapdf_pdf, double Q0, double Qf)
		: _pdf{std::move(lhapdf_pdf)}
	{
		_Q0 = Q0;

		// up and down (1, 2) are going to be considered massless always (for now)
		for (uint i : {3, 4, 5, 6})
			_masses[i] = _pdf->quarkMass(i);

		// set the first mass not less than Qf to Qf,
		// and zero everything above it
		auto it = std::lower_bound(_masses.begin(), _masses.end(), Qf);
		if (it == _masses.end())
			log(LOG_ERROR, "LHAPDF", "std::lower_bound failed");
		*it =  Qf;
		while (++it != _masses.end())
			*it = 0.0;

		// from reverse, set first mass less than Q0 to Q0,
		// and everything after it to zero
		it = std::upper_bound(_masses.begin(), _masses.end(), Q0);
	    if (it == _masses.end()) {
			log(LOG_ERROR, "LHAPDF", "std::upper_bound failed");
		} else if (it == _masses.begin()) {
			log(LOG_ERROR, "LHAPDF", "invalid provided mass array. input energy too low");
		}

		*(--it) = Q0;
		_nfi = std::distance(_masses.begin(), it);
		while (--it != _masses.begin())
			*it = 0.0;
		
		_alpha0 = _pdf->alphasQ(Q0);

		std::ostringstream ss{};
		std::copy(_masses.begin(), _masses.end(), std::ostream_iterator<double>(ss, ", "));
		log(LOG_DEBUG, "LHAPDF", "Using the following masses: {}", ss.str());
		log(LOG_DEBUG, "LHAPDF", "Using nfi = {}", _nfi);
		log(LOG_DEBUG, "LHAPDF", "Using alpha0 = {}", _alpha0);
		
	}

	void LHAPDFDistribution::fillSingletCoeffs(
		coefficient_accessor_type const& accessor,
		std::vector<double> const& grid_points) const
	{
		for (uint k=0; k<grid_points.size()-1; k++) {
			double x = grid_points[k];
			accessor(0, k) = xg(x);
			accessor(1, k) = xqplus(x);
		}
	}
	
	void LHAPDFDistribution::fillNonSingletCoeffs(
		coefficient_accessor_type const& accessor,
		std::vector<double> const& grid_points) const
	{
		for (uint k=0; k<grid_points.size()-1; k++) {
			double x = grid_points[k];
			accessor(1, k) = xu(x);  // u
			accessor(2, k) = xd(x);  // d
			accessor(3, k) = xs(x);  // s
			accessor(7, k) = xub(x); // ub
			accessor(8, k) = xdb(x); // db
			accessor(9, k) = xs(x);  // sb ( = s)
		}
	}
}
