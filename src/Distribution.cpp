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


	void LesHouchesQED::fillCoeffs(
		accessor_type const& s_accessor,
		accessor_type const& ns_accessor,
		std::vector<double> const& grid_points) const
	{
		for (uint k=0; k<grid_points.size()-1; k++) {
			double x = grid_points[k];

			s_accessor(static_cast<uint>(QEDPartonIndices::SIGMAUD), k) = xsigmaud(x);
			s_accessor(static_cast<uint>(QEDPartonIndices::SIGMA), k) = xsigma(x); 
			s_accessor(static_cast<uint>(QEDPartonIndices::G), k) = xg(x);
			s_accessor(static_cast<uint>(QEDPartonIndices::PHOTON), k) = xgamma(x);
			s_accessor(static_cast<uint>(QEDPartonIndices::SIGMAL), k) = xsigmal(x);
			
			ns_accessor(static_cast<uint>(QEDPartonIndices::U), k) = xu(x);  
			ns_accessor(static_cast<uint>(QEDPartonIndices::D), k) = xd(x);  
			ns_accessor(static_cast<uint>(QEDPartonIndices::S), k) = xs(x);  
			ns_accessor(static_cast<uint>(QEDPartonIndices::UB), k) = xub(x);
			ns_accessor(static_cast<uint>(QEDPartonIndices::DB), k) = xdb(x);
			ns_accessor(static_cast<uint>(QEDPartonIndices::SB), k) = xs(x);

			ns_accessor(static_cast<uint>(QEDPartonIndices::UV), k) = xuv(x);  
			ns_accessor(static_cast<uint>(QEDPartonIndices::DV), k) = xdv(x);

			ns_accessor(static_cast<uint>(QEDPartonIndices::SIGMAUC), k) = xsigmauc(x);
			ns_accessor(static_cast<uint>(QEDPartonIndices::SIGMADS), k) = xsigmads(x);
			ns_accessor(static_cast<uint>(QEDPartonIndices::SIGMASB), k) = xsigmasb(x);
			
			ns_accessor(static_cast<uint>(QEDPartonIndices::E), k) = xe(x); 
			ns_accessor(static_cast<uint>(QEDPartonIndices::MU), k) = xmu(x); 
			ns_accessor(static_cast<uint>(QEDPartonIndices::TAU), k) = xtau(x); 
			ns_accessor(static_cast<uint>(QEDPartonIndices::EB), k) = xeb(x);
			ns_accessor(static_cast<uint>(QEDPartonIndices::MUB), k) = xmub(x);
			ns_accessor(static_cast<uint>(QEDPartonIndices::TAUB), k) = xtaub(x); 
		}
	}

	void LesHouchesQED::setup(double q0, double qf)
	{
		_Q0 = q0;
		_Qf = qf;
		_alpha0 = 0.35;
		_alphaqed0 = 1.0/133.4;
		_nfi = 4;
		_nff = 4;
	}


} // namespace Candia2
