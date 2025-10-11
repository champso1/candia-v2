#include "Candia-v2/SplittingFn.hpp"
#include "Candia-v2/SpecialFuncs.hpp"

namespace Candia2
{
	namespace splitfunc_fortran
    {

		double P2nsm::Regular(const double _x) const
		{
			// avoid const shenanigans
			// `T const*` cannot be converted to `T *`
			
			double x = _x;
			int nf = static_cast<int>(_nf);
			
			return p2nsma_(&x, &nf)/8.0;
		}
		double P2nsm::Plus(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2nsb_(&x, &nf)/8.0;
		}
		double P2nsm::Delta(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2nsmc_(&x, &nf)/8.0;
		}


		double P2nsp::Regular(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2nspa_(&x, &nf)/8.0;
		}
		double P2nsp::Plus(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2nsb_(&x, &nf)/8.0;
		}
		double P2nsp::Delta(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2nspc_(&x, &nf)/8.0;
		}

		double P2nsv::Regular(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return (p2nsma_(&x, &nf) + p2nssa_(&x, &nf))/8.0;
		}
		double P2nsv::Plus(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2nsb_(&x, &nf)/8.0;
		}
		double P2nsv::Delta(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2nsmc_(&x, &nf)/8.0;
		}

		double P2qq::Regular(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return (p2nspa_(&x, &nf) + p2psa_(&x, &nf))/8.0;
		}
		double P2qq::Plus(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2nsb_(&x, &nf)/8.0;
		}
		double P2qq::Delta(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2nspc_(&x, &nf)/8.0;
		}

		double P2qg::Regular(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2qga_(&x, &nf)/8.0;
		}

		double P2gq::Regular(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2gqa_(&x, &nf)/8.0;
		}

		double P2gg::Regular(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2gga_(&x, &nf)/8.0;
		}
		double P2gg::Plus(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2ggb_(&x, &nf)/8.0;
		}
		double P2gg::Delta(const double _x) const
		{
			double x = _x;
			int nf = static_cast<int>(_nf);

			return p2ggc_(&x, &nf)/8.0;
		}
		
	} // namespace hplog
} // namespace Candia2
