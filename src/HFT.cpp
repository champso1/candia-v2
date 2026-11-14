// HFT.cpp

#include "Candia-v2/ArrayGrid.hpp"
#include "Candia-v2/Candia.hpp"
#include "Parallel.hpp"

#include <algorithm>
#include <execution>

namespace Candia2
{
    void DGLAPSolver::heavyFlavorTreatment()
    {
		if (options.disable_heavy_flavor_matching) {
			log(LOG_WARNING, "HFT", "Matching has been disabled.");
			return;
		}
		
		log(LOG_DEBUG, "HFT", "Treating heavy flavors: {}th quark mass threshold (mass {})", _nf+1, _alpha_s.masses(_nf+1));
		// Copy of pre-threshold distributions
		// the nf+1 dists are defined in terms of the nf dists,
		// so we need this copy since we would otherwise be overwriting
		// a given dist as we do the calculations
		// we just store them in the s=1 array, and modify the s=0 array
		// (which are the initial conditions for the next set of iterations
		// at the next nf
		log(LOG_DEBUG, "HFT", "Creating copy of pre-threshold distributions... ");
        
        std::vector<ArrayGrid> arr(13, ArrayGrid(_grid.size()));
		std::vector<ArrayGrid> arr_singlet(2, ArrayGrid(_grid.size()));

		for (uint j=0; j<=1; ++j)
			std::ranges::copy(_S(0,j,0), arr_singlet[j].begin());

		auto arr_accessor = [&](uint j) -> ArrayGridView {
			if (_order == 2) {
				return _evol_type == EvolType::Truncated ? _S_NS(0,j,0) : _C(j,0,0,0);
			} else if (_order == 3) {
				return _evol_type == EvolType::Truncated ? _S_NS(0,j,0) : _D(j,0,0,0,0);
			}
			throw std::runtime_error("unreachable");
		};
		
		for (uint i=1; i<=_nf; i++) {
			for (uint j=i; j<=i+6; j+=6)
				std::ranges::copy(arr_accessor(j), arr[j].begin());
		}

		double as = _alpha_s.post(_nf+1);
		log(LOG_DEBUG, "HFT", "Value of alpha_s post threshold: {}", as);

		// precompute the convolutions with the OMEs used below
		std::vector<ExprName> omes{ExprName::A2ns, ExprName::A2gq, ExprName::A2gg, ExprName::A2hq, ExprName::A2hg};
		if (_order == 3 && !options.use_nnlo_matching_conditions_at_n3lo) {
			omes.insert(omes.end(), {ExprName::A3nsm, ExprName::A3nsp, ExprName::A3gq, ExprName::A3gg,
				ExprName::A3hq, ExprName::A3hg, ExprName::A3psqq, ExprName::A3sqg});
			if (options.use_n3lo_heavyquark_asymmetry)
				omes.push_back(ExprName::A3PSshq);
		}
		for (ExprName name : omes)
			buildKernel(name);

		auto grid_idxs = std::ranges::views::iota(uint{0}, _grid.size());
		
		if (_order == 2) {
			parallelFor(grid_idxs, [&](uint k) {
				// q
				for (uint j=1; j<=_nf; j++)
					HFT_NNLO1(arr[j], k, arr_accessor(j));
				// qbar
				for (uint j=1+6; j<=_nf+6; j++)
					HFT_NNLO1(arr[j], k, arr_accessor(j));

				HFT_NNLO2(arr_singlet[0], arr_singlet[1], k); // gluon
				HFT_NNLO3(arr_singlet[0], arr_singlet[1], k, arr_accessor(_nf+1), arr_accessor(_nf+1+6)); // heavy flavor
			});
		} else if (_order == 3) {
			if (!options.use_nnlo_matching_conditions_at_n3lo) {
				log(LOG_DEBUG, "HFT", "Performing N3LO matching at N3LO");
				ArrayGrid qminus(_grid.size());
				for (uint k=0; k<_grid.size(); ++k) {
					qminus[k] = 0.0;
					for (uint j=1; j<=_nf; ++j) {
						qminus[k] += arr[j][k] - arr[j+6][k];
					}
				}

				auto& a3psqq = getKernel(ExprName::A3psqq);
				auto& a3sqg  = getKernel(ExprName::A3sqg);
				const double fac_n3lo = as*as*as/(64.0*PI_3);
				parallelFor(grid_idxs, [&](uint k) {
					const double convSPa = a3psqq.apply(arr_singlet[1], k);
					const double convSPb = a3sqg.apply(arr_singlet[0], k);
					const double SP = fac_n3lo*(convSPa + convSPb)/static_cast<double>(_nf);

					// q
					for (uint j=1; j<=_nf; j++)
						HFT_N3LO1(arr[j], arr[j+6], k, SP, arr_accessor(j));
					// qbar
					for (uint j=1+6; j<=_nf+6; j++)
						HFT_N3LO2(arr[j-6], arr[j], k, SP, arr_accessor(j));

					HFT_N3LO3(arr_singlet[0], arr_singlet[1], k); // gluon
					HFT_N3LO4(arr_singlet[0], arr_singlet[1], qminus, k, arr_accessor(_nf+1), arr_accessor(_nf+1+6)); // heavy flavor
				});
			} else {
				log(LOG_DEBUG, "HFT", "Performing NNLO matching at N3LO");
				parallelFor(grid_idxs, [&](uint k) {
					// q
					for (uint j=1; j<=_nf; j++)
						HFT_NNLO1(arr[j], k, arr_accessor(j));
					// qbar
					for (uint j=1+6; j<=_nf+6; j++)
						HFT_NNLO1(arr[j], k, arr_accessor(j));

					HFT_NNLO2(arr_singlet[0], arr_singlet[1], k); // gluon
					HFT_NNLO3(arr_singlet[0], arr_singlet[1], k, arr_accessor(_nf+1), arr_accessor(_nf+1+6)); // heavy flavor
				});
			}
		}

		// the OMEs change with the number of flavours, so these kernels are not reused
		for (ExprName name : omes)
			_kernels[static_cast<uint>(name)] = ConvolutionKernel{};
    }

	void DGLAPSolver::HFT_NNLO1(ArrayGridView c, uint k, ArrayGridView q)
    {
		auto& a2ns = getKernel(ExprName::A2ns);
        double const as = _alpha_s.post(_nf+1);
        double const conv = a2ns.apply(c, k);
        
	    q[k] += std::pow(as/(4.0*PI), 2) * conv;
    }

    void DGLAPSolver::HFT_NNLO2(ArrayGridView g, ArrayGridView qp, uint k)
    {
		auto& a2gq = getKernel(ExprName::A2gq);
		auto& a2gg = getKernel(ExprName::A2gg);
        double const as = _alpha_s.post(_nf+1);
        double const conv1 = a2gq.apply(qp, k);
        double const conv2 = a2gg.apply(g, k);

		_S(0,0,0,k) += std::pow(as/(4.0*PI), 2) * (conv1 + conv2);
    }

    void DGLAPSolver::HFT_NNLO3(ArrayGridView g, ArrayGridView qp, uint k, ArrayGridView qh, ArrayGridView qhb)
    {
		auto& a2hq = getKernel(ExprName::A2hq);
		auto& a2hg = getKernel(ExprName::A2hg);
        double const as = _alpha_s.post(_nf+1);
        double const conv1 = a2hq.apply(qp, k);
        double const conv2 = a2hg.apply(g, k);
		double const fac = 0.5*std::pow(as/(4.0*PI), 2) * (conv1 + conv2);

		qh[k] = fac;
        qhb[k] = fac;
    }

    // q
	void DGLAPSolver::HFT_N3LO1(ArrayGridView q, ArrayGridView qb, uint k, double SP, ArrayGridView qh)
	{
	    const double as = _alpha_s.post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

		double const conv1a = getKernel(ExprName::A2ns).apply(q, k);
		double const conv1b = getKernel(ExprName::A3nsp).apply(q, k);

		double const conv2a = getKernel(ExprName::A2ns).apply(qb, k);
		double const conv2b = getKernel(ExprName::A3nsp).apply(qb, k);

		double const conv3a = conv1a;
		double const conv3b = getKernel(ExprName::A3nsm).apply(q, k);

		double const conv4a = conv2a;
		double const conv4b = getKernel(ExprName::A3nsm).apply(qb, k);

		qh[k] += 0.5*(
			((fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b)) +
			((fac_nnlo*conv3a + fac_n3lo*conv3b) - (fac_nnlo*conv4a + fac_n3lo*conv4b)) +
			SP);
	}

    // qbar
	void DGLAPSolver::HFT_N3LO2(ArrayGridView q, ArrayGridView qb, uint k, double SP, ArrayGridView qhb)
	{
		const double as = _alpha_s.post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

		double const conv1a = getKernel(ExprName::A2ns).apply(q, k);
		double const conv1b = getKernel(ExprName::A3nsp).apply(q, k);

		double const conv2a = getKernel(ExprName::A2ns).apply(qb, k);
		double const conv2b = getKernel(ExprName::A3nsp).apply(qb, k);

		double const conv3a = conv1a;
		double const conv3b = getKernel(ExprName::A3nsm).apply(q, k);

		double const conv4a = conv2a;
		double const conv4b = getKernel(ExprName::A3nsm).apply(qb, k);

		qhb[k] += 0.5*(
			((fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b)) -
			((fac_nnlo*conv3a + fac_n3lo*conv3b) - (fac_nnlo*conv4a + fac_n3lo*conv4b)) +
			SP);
	}

	// gluon (index 0 in S array)
	void DGLAPSolver::HFT_N3LO3(ArrayGridView g, ArrayGridView qp, uint k)
	{
		const double as = _alpha_s.post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);
	    
		const double conv1a = getKernel(ExprName::A2gq).apply(qp, k);
		const double conv1b = getKernel(ExprName::A3gq).apply(qp, k);
		const double conv2a = getKernel(ExprName::A2gg).apply(g, k);
		const double conv2b = getKernel(ExprName::A3gg).apply(g, k);
		
		_S(0,0,0,k) += (fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b);
	}

	// heavy flavor
	void DGLAPSolver::HFT_N3LO4(ArrayGridView g, ArrayGridView qp, ArrayGridView qminus, uint k, ArrayGridView qh, ArrayGridView qhb)
	{
		const double as = _alpha_s.post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

	    const double conv1a = getKernel(ExprName::A2hq).apply(qp, k);
		const double conv1b = getKernel(ExprName::A3hq).apply(qp, k);
		const double conv2a = getKernel(ExprName::A2hg).apply(g, k);
		const double conv2b = getKernel(ExprName::A3hg).apply(g, k);
		const double conv3  =
			options.use_n3lo_heavyquark_asymmetry ?
			getKernel(ExprName::A3PSshq).apply(qminus, k)
			: 0.0;

		const double res = (fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b);
        const double res_qh = 0.5*(res + fac_n3lo*conv3);
		const double res_qb = 0.5*(res - fac_n3lo*conv3);
		qh[k] = res_qh;
        qhb[k] = res_qb;
	}
	
} // namespace Candia2
