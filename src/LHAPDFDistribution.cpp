// LHAPDFDistribution.cpp

#include "Candia-v2/LHAPDFDistribution.hpp"

namespace Candia2
{	
	void LHAPDFDistribution::fillCoeffs(
		accessor_type const& s_accessor,
		accessor_type const& ns_accessor,
		std::vector<value_type> const& grid_points) const
	{
		for (uint k=0; k<grid_points.size()-1; k++) {
			double x = grid_points[k];

			s_accessor(0, k) = xg(x);
			s_accessor(1, k) = xqplus(x);
			
			ns_accessor(1, k) = xu(x);  // u
			ns_accessor(2, k) = xd(x);  // d
			ns_accessor(3, k) = xs(x);  // s
			ns_accessor(7, k) = xub(x); // ub
			ns_accessor(8, k) = xdb(x); // db
			ns_accessor(9, k) = xs(x);  // sb ( = s)
		}
	}

	void LHAPDFDistribution::setup(double Q0, double Qf)
	{
		if (Q0 < 1.0)
			log(LOG_ERROR, "LHAPDFDistribution", "Q0 cannot be less than 1 GeV");
		
		_Q0 = Q0;
		_Qf = Qf;

		// up and down (1, 2) are going to be considered massless
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
			log(LOG_ERROR, "LHAPDFDistribution", "error finding nfi (_nfi={})", _nfi);

		for (_nff=6; _nff>=_nfi; --_nff) {
			if (Qf > _masses[_nff])
				break;
		}
		// we don't need the equivalent second if statement here like with nfi
		// because again, say, if charm mass (nf=4) = 1.3 and bottom mass (nf=5) = 4.5,
		// if we have Qf=2, or anything up to and including 4.5,
		// then we are good staying with _nff = 4
		if (_nff < 1 || _nff > 6)
			log(LOG_ERROR, "LHAPDFDistribution", "error finding nf (_nff={})", _nff);

		// TODO: should we determine it nicely from the reference scale?
	    _alpha0 = _pdf->alphasQ(Q0);

		log(LOG_DEBUG, "LHAPDFDistribution", "Using the following masses: {}", vec_to_str(_masses));
		log(LOG_DEBUG, "LHAPDFDistribution", "Using nfi = {}", _nfi);
		log(LOG_DEBUG, "LHAPDFDistribution", "Using nff = {}", _nff);
		log(LOG_DEBUG, "LHAPDFDistribution", "Using alpha0 = {} at q0 = {}", _alpha0, _Q0);
	}
	
} // namespace Candia2
