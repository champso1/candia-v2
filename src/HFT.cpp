#include "Candia-v2/Candia.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"

namespace Candia2
{
    void DGLAPSolver::heavyFlavorTreatment()
    {
		if (options.disable_heavy_flavor_matching) {
			log(LOG_WARNING, "HFT", "Matching has been disabled.");
			return;
		}
		
		log(LOG_INFO, "HFT", "Treating heavy flavors: {}th quark mass threshold (mass {})", _nf+1, _alpha_s.masses(_nf+1));
		OpMatElem::update(-_log_mur2_muf2, _nf);

		// Copy of pre-threshold distributions
		// the nf+1 dists are defined in terms of the nf dists,
		// so we need this copy since we would otherwise be overwriting
		// a given dist as we do the calculations
		// we just store them in the s=1 array, and modify the s=0 array
		// (which are the initial conditions for the next set of iterations
		// at the next nf
		log(LOG_INFO, "HFT", "Creating copy of pre-threshold distributions... ");
        
        std::vector<ArrayGrid> arr(13, ArrayGrid(_grid.size()));
		std::vector<ArrayGrid> arr_singlet(2, ArrayGrid(_grid.size()));

		for (uint j=0; j<=1; ++j)
			arr_singlet[j] = _S[0][j][0];

		auto arr_accessor = [&](uint j) -> ArrayGrid& {
			if (_order == 2) {
				return options.use_truncated_nonsinglet_sol ? _S_NS[0][j][0] : _C[j][0][0][0];
			} else if (_order == 3) {
				return options.use_truncated_nonsinglet_sol ? _S_NS[0][j][0] : _D[j][0][0][0][0];
			}
			throw "unreachable";
		};
		
		if (_order == 2) {
			for (uint i=1; i<=_nf; i++) {
				for (uint j=i; j<=i+6; j+=6)
					arr[j] = arr_accessor(j);
			}
		} else if (_order == 3) {
            for (uint i=1; i<=_nf; i++) {
				for (uint j=i; j<=i+6; j+=6)
					arr[j] = arr_accessor(j);
			}
        }

		double as = _alpha_s.post(_nf+1);
		log(LOG_INFO, "HFT", "Value of alpha_s post threshold: {}", as);

		if (_order == 2) {
			for (uint k=0; k<_grid.size()-1;k++) {
				// q
				for (uint j=1; j<=_nf; j++)
					HFT_NNLO1(arr[j], k, arr_accessor(j));
				// qbar
				for (uint j=1+6; j<=_nf+6; j++)
					HFT_NNLO1(arr[j], k, arr_accessor(j));

				HFT_NNLO2(arr_singlet[0], arr_singlet[1], k); // gluon
				HFT_NNLO3(arr_singlet[0], arr_singlet[1], k, arr_accessor(_nf+1), arr_accessor(_nf+1+6)); // heavy flavor
			}
		} else if (_order == 3) {
			if (!options.use_nnlo_matching_conditions_at_n3lo) {
				log(LOG_INFO, "HFT", "Performing N3LO matching at N3LO");
				ArrayGrid qminus(_grid.size());
				for (uint k=0; k<_grid.size(); ++k) {
					qminus[k] = 0.0;
					for (uint j=1; j<=_nf; ++j) {
						qminus[k] += arr[j][k] - arr[j+6][k];
					}
				}

				static auto& a3psqq = getExpression("A3psqq");
				static auto& a3sqg  = getExpression("A3sqg");
				for (uint k=0; k<_grid.size()-1;k++) {
					const double fac_n3lo = as*as*as/(64.0*PI_3);
					const double convSPa = _grid.convolution(arr_singlet[1], a3psqq, k);
					const double convSPb = _grid.convolution(arr_singlet[0], a3sqg, k);
					const double SP = fac_n3lo*(convSPa + convSPb)/static_cast<double>(_nf);

					// q
					for (uint j=1; j<=_nf; j++)
						HFT_N3LO1(arr[j], arr[j+6], j, k, SP, arr_accessor(j));
					// qbar
					for (uint j=1+6; j<=_nf+6; j++)
						HFT_N3LO2(arr[j-6], arr[j], j, k, SP, arr_accessor(j));

					HFT_N3LO3(arr_singlet[0], arr_singlet[1], k); // gluon
					HFT_N3LO4(arr_singlet[0], arr_singlet[1], qminus, k, arr_accessor(_nf+1), arr_accessor(_nf+1+6)); // heavy flavor
				}
			} else {
				log(LOG_INFO, "HFT", "Performing NNLO matching at N3LO");
				for (uint k=0; k<_grid.size()-1;k++) {
					// q
					for (uint j=1; j<=_nf; j++)
						HFT_NNLO1(arr[j], k, arr_accessor(j));
					// qbar
					for (uint j=1+6; j<=_nf+6; j++)
						HFT_NNLO1(arr[j], k, arr_accessor(j));

					HFT_NNLO2(arr_singlet[0], arr_singlet[1], k); // gluon
					HFT_NNLO3(arr_singlet[0], arr_singlet[1], k, arr_accessor(_nf+1), arr_accessor(_nf+1+6)); // heavy flavor
				}
			}
		}
    }
    void DGLAPSolver::HFT_NNLO1(ArrayGrid& c, uint k, ArrayGrid& q)
    {
		static auto& a2ns = getExpression("A2ns");
        double const as = _alpha_s.post(_nf+1);
        double const conv = _grid.convolution(c, a2ns, k);
        
	    q[k] += std::pow(as/(4.0*PI), 2) * conv;
    }
    void DGLAPSolver::HFT_NNLO2(ArrayGrid& g, ArrayGrid& qp, uint k)
    {
		static auto& a2gq = getExpression("A2gq");
		static auto& a2gg = getExpression("A2gg");
        double const as = _alpha_s.post(_nf+1);
        double const conv1 = _grid.convolution(qp, a2gq, k);
        double const conv2 = _grid.convolution(g, a2gg, k);

		_S[0][0][0][k] += std::pow(as/(4.0*PI), 2) * (conv1 + conv2);
    }
    void DGLAPSolver::HFT_NNLO3(ArrayGrid& g, ArrayGrid& qp, uint k, ArrayGrid& qh, ArrayGrid& qhb)
    {
		static auto& a2hq = getExpression("A2hq");
		static auto& a2hg = getExpression("A2hg");
        double const as = _alpha_s.post(_nf+1);
        double const conv1 = _grid.convolution(qp, a2hq, k);
        double const conv2 = _grid.convolution(g, a2hg, k);
		double const fac = 0.5*std::pow(as/(4.0*PI), 2) * (conv1 + conv2);

		qh[k] = fac;
        qhb[k] = fac;
    }

    // q
	void DGLAPSolver::HFT_N3LO1(ArrayGrid& q, ArrayGrid& qb, uint j, uint k, double SP, ArrayGrid& qh)
	{
	    const double as = _alpha_s.post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

		double const conv1a = _grid.convolution(q, getExpression("A2ns"), k);
		double const conv1b = _grid.convolution(q, getExpression("A3nsp"), k);

		double const conv2a = _grid.convolution(qb, getExpression("A2ns"), k);
		double const conv2b = _grid.convolution(qb, getExpression("A3nsp"), k);

		double const conv3a = conv1a;
		double const conv3b = _grid.convolution(q, getExpression("A3nsm"), k);

		double const conv4a = conv2a;
		double const conv4b = _grid.convolution(qb, getExpression("A3nsm"), k);

		qh[k] += 0.5*(
			((fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b)) +
			((fac_nnlo*conv3a + fac_n3lo*conv3b) - (fac_nnlo*conv4a + fac_n3lo*conv4b)) +
			SP);
	}

    // qbar
	void DGLAPSolver::HFT_N3LO2(ArrayGrid& q, ArrayGrid& qb, uint j, uint k, double SP, ArrayGrid& qhb)
	{
		const double as = _alpha_s.post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

		double const conv1a = _grid.convolution(q, getExpression("A2ns"), k);
		double const conv1b = _grid.convolution(q, getExpression("A3nsp"), k);

		double const conv2a = _grid.convolution(qb, getExpression("A2ns"), k);
		double const conv2b = _grid.convolution(qb, getExpression("A3nsp"), k);

		double const conv3a = conv1a;
		double const conv3b = _grid.convolution(q, getExpression("A3nsm"), k);

		double const conv4a = conv2a;
		double const conv4b = _grid.convolution(qb, getExpression("A3nsm"), k);

		qhb[k] += 0.5*(
			((fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b)) -
			((fac_nnlo*conv3a + fac_n3lo*conv3b) - (fac_nnlo*conv4a + fac_n3lo*conv4b)) +
			SP);
	}

	// gluon (index 0 in S array)
	void DGLAPSolver::HFT_N3LO3(ArrayGrid& g, ArrayGrid& qp, uint k)
	{
		const double as = _alpha_s.post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);
	    
		const double conv1a = _grid.convolution(qp, getExpression("A2gq"), k);
		const double conv1b = _grid.convolution(qp, getExpression("A3gq"), k);
		const double conv2a = _grid.convolution(g, getExpression("A2gg"), k);
		const double conv2b = _grid.convolution(g, getExpression("A3gg"), k);
		
		_S[0][0][0][k] += (fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b);
	}

	// heavy flavor
	void DGLAPSolver::HFT_N3LO4(ArrayGrid& g, ArrayGrid& qp, ArrayGrid& qminus, uint k, ArrayGrid& qh, ArrayGrid& qhb)
	{
		const double as = _alpha_s.post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

	    const double conv1a = _grid.convolution(qp, getExpression("A2hq"), k);
		const double conv1b = _grid.convolution(qp, getExpression("A3hq"), k);
		const double conv2a = _grid.convolution(g, getExpression("A2hg"), k);
		const double conv2b = _grid.convolution(g, getExpression("A3hg"), k);
		const double conv3  =
			options.use_n3lo_heavyquark_asymmetry ?
			_grid.convolution(qminus, getExpression("AQqPSs3"), k)
			: 0.0;

		const double res = (fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b);
        const double res_qh = 0.5*(res + fac_n3lo*conv3);
		const double res_qb = 0.5*(res - fac_n3lo*conv3);
		qh[k] = res_qh;
        qhb[k] = res_qb;
	}
} // namespace Candia2
