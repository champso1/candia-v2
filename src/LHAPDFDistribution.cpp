// LHAPDFDistribution.cpp

#include "Candia-v2/LHAPDFDistribution.hpp"

namespace Candia2
{
	LHAPDFDistribution::LHAPDFDistribution(lhapdf_pdf_ptr_type lhapdf_pdf, value_type Q0, value_type Qf)
		: _pdf{std::move(lhapdf_pdf)}
	{
		// 
		_Q0 = Q0;
		_Qf = Qf;

		// TODO: not sure how I feel interpolating this value of alphas
		// then proceeding to re-evolve it
		// we should probably either determine this in some good way from the
		// PDFs reference scale, or let alphas itself be grabbed from the PDF
		// during the evolution
		_alpha0 = _pdf->alphasQ(Q0);
		
		// up and down (1, 2) are going to be considered massless always (for now)
		for (uint i : {3, 4, 5, 6})
			_masses[i] = _pdf->quarkMass(i);

		for (_nfi=1; _nfi<=_masses.size(); ++_nfi) {
			if (Q0 <= _masses[_nfi])
				break;
		}
		// this is required because if, say, charm mass (nf=4) is 1.3
		// and q0 = 1.0, the above loop sets _nfi=4, but in reality
		// it should be nf=3 since q0<1.3
		if (Q0 < _masses[_nfi])
			_nfi--;
		if (_nfi < 1 || _nfi > 6)
			log(LOG_ERROR, "LHAPDF", "error finding nfi (_nfi={})", _nfi);

		for (_nff=6; _nff>=_nfi; --_nff) {
			if (Qf > _masses[_nff])
				break;
		}
		// we don't need the equivalent second if statement here like with nfi
		// because again, say, if charm mass (nf=4) = 1.3 and bottom mass (nf=5) = 4.5,
		// if we have Qf=2, or anything up to and including 4.5,
		// then we are good staying with _nff = 4
		if (_nff < 1 || _nff > 6)
			log(LOG_ERROR, "LHAPDF", "error finding nf (_nff={})", _nff);

		log(LOG_DEBUG, "LHAPDF", "Using the following masses: {}", vec_to_str(_masses));
		log(LOG_DEBUG, "LHAPDF", "Using nfi = {}", _nfi);
		log(LOG_DEBUG, "LHAPDF", "Using nff = {}", _nff);
		log(LOG_DEBUG, "LHAPDF", "Using alpha0 = {}", _alpha0);
		
	}

	void LHAPDFDistribution::fillSingletCoeffs(
		accessor_type const& accessor,
		std::vector<value_type> const& grid_points) const
	{
		for (uint k=0; k<grid_points.size()-1; k++) {
			double x = grid_points[k];
			accessor(0, k) = xg(x);
			accessor(1, k) = xqplus(x);
		}
	}
	
	void LHAPDFDistribution::fillNonSingletCoeffs(
		accessor_type const& accessor,
		std::vector<value_type> const& grid_points) const
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
} // namespace Candia2
