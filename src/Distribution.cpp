// Distribution.cpp

#include "Candia-v2/Distribution.hpp"

namespace Candia2
{
	LesHouchesDistribution::LesHouchesDistribution(double qf)
	{
		_Q0 = std::numbers::sqrt2;
		_Qf = qf;
		_alpha0 = 0.35;
		_nfi = 3;
		_masses = _leshouche_masses;

		for (_nff=6; _nff>=_nfi; --_nff) {
			if (qf > _masses[_nff])
				break;
		}
		if (_nff < 1 || _nff > 6)
			log(LOG_ERROR, "LHAPDF", "error finding nf (_nff={})", _nff);
	}
	
	void LesHouchesDistribution::fillSingletCoeffs(
		accessor_type const& accessor,
		std::vector<value_type> const& grid_points) const
	{
		for (uint k=0; k<grid_points.size()-1; k++) {
			double x = grid_points[k];
			accessor(0, k) = xg(x);
			accessor(1, k) = xqplus(x);
		}
	}
	
	void LesHouchesDistribution::fillNonSingletCoeffs(
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
