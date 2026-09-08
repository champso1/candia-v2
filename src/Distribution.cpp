// Distribution.cpp

#include "Candia-v2/Distribution.hpp"

namespace Candia2
{
	void LesHouchesDistribution::fillCoeffs(
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

	void LesHouchesDistribution::setup(double q0, double qf)
	{
		_Q0 = q0;
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


	void QEDDistribution::fillCoeffs(
		accessor_type const& s_accessor,
		accessor_type const& ns_accessor,
		std::vector<double> const& grid_points) const
	{
		for (uint k=0; k<grid_points.size()-1; k++) {
			double x = grid_points[k];

			// this is a clusterfuck
			double u = xuv(x) + 0.9*xub(x);
			double d = xdv(x) + 0.9*xdb(x);
			double s = xs(x);
			double ub = 0.9*xub(x);
			double db = 0.9*xdb(x);
			double sb = xsb(x);
			double g = xg(x)*0.99;

			double sigma =
				u + ub +
				d + db +
				s + sb;
			double ep = xdb(x)/10.0;
			double mup = xub(x)/10.0;
			double taup = xdb(x)/10.0;
			double sigmal = ep + mup + taup;
			double gamma = xg(x)/100.0 + xub(x)/10.0;
			double deltaUD =
				u + ub +
				- d - db
				- s - sb;
			double deltads =
				d + db
				- s - sb;
			double deltal2 = ep - mup;
			double deltal3 = ep + mup - 2.0*taup;

			s_accessor(static_cast<uint>(QEDPartonIndices::DELTAUD), k) = deltaUD;
			s_accessor(static_cast<uint>(QEDPartonIndices::SIGMA), k) = sigma; 
			s_accessor(static_cast<uint>(QEDPartonIndices::G), k) = g;
			s_accessor(static_cast<uint>(QEDPartonIndices::PHOTON), k) = gamma;
			s_accessor(static_cast<uint>(QEDPartonIndices::SIGMAL), k) = sigmal;

			ns_accessor(static_cast<uint>(QEDPartonIndices::UV), k) = xuv(x);  
			ns_accessor(static_cast<uint>(QEDPartonIndices::DV), k) = xdv(x);
			ns_accessor(static_cast<uint>(QEDPartonIndices::DELTADS), k) = deltads;
			ns_accessor(static_cast<uint>(QEDPartonIndices::DELTAL2), k) = deltal2;
			ns_accessor(static_cast<uint>(QEDPartonIndices::DELTAL3), k) = deltal3;
		}
	}

	void QEDDistribution::setup(double q0, double qf)
	{
		_Q0 = q0;
		_Qf = qf;
		_alpha0 = 0.35;
		_alphaqed0 = 0.008539696327675218;
		_nfi = 4;
		_nff = 4;
	}


} // namespace Candia2
