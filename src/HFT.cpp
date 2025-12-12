#include "Candia-v2/Candia.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"

#include <print>

namespace Candia2
{
    void DGLAPSolver::heavyFlavorTreatment()
    {
		std::println("[DGLAP] Treating heavy flavors: {}th quark mass threshold (mass {})", _nf+1, _alpha_s.masses(_nf+1));
		OpMatElem::update(-_log_mur2_muf2, _nf);

		// Copy of pre-threshold distributions
		// the nf+1 dists are defined in terms of the nf dists,
		// so we need this copy since we would otherwise be overwriting
		// a given dist as we do the calculations
		// we just store them in the s=1 array, and modify the s=0 array
		// (which are the initial conditions for the next set of iterations
		// at the next nf
		std::print("[DGLAP] Creating copy of pre-threshold distributions... ");
        
        std::vector<ArrayGrid> arr(13, ArrayGrid{_grid});
        std::vector<ArrayGrid> arr_singlet(2, ArrayGrid{_grid});

		for (uint j=0; j<=1; ++j)
			arr_singlet[j] = _S2[0][j][0];
		
		if (_order == 2) 
		{
			for (uint i=1; i<=_nf; i++)
			{
				for (uint j=i; j<=i+6; j+=6)
					arr[j] = _C2[j][0][0][0];
			}
		} else if (_order == 3) {
            for (uint i=1; i<=_nf; i++)
			{
				for (uint j=i; j<=i+6; j+=6)
					arr[j] = _D2[j][0][0][0][0];
			}
        }
		std::println("Done.");

		double as = _alpha_s.post(_nf+1);
		std::println("[DGLAP] Value of alpha_s post threshold: {}", as);

		switch(_order)
		{
			case 2:
			{
				for (uint k=0; k<_grid.size()-1;k++)
				{
					// q
					for (uint j=1; j<=_nf; j++)
						HFT_NNLO1(arr[j], j, k);
					// qbar
					for (uint j=1+6; j<=_nf+6; j++)
						HFT_NNLO1(arr[j], j, k);

					HFT_NNLO2(arr_singlet[1], arr_singlet[0], k); // gluon
					HFT_NNLO3(arr_singlet[1], arr_singlet[0], k); // heavy flavor
				}
			} break;
            case 3:
            {
                for (uint k=0; k<_grid.size()-1;k++)
				{
                    const double fac_n3lo = as*as*as/(64.0*PI_3);
					const double convSPa = getSplitFunc("A3psqq").convolution(arr_singlet[1], k);
					const double convSPb = getSplitFunc("A3sqg").convolution(arr_singlet[0], k);
					const double SP = fac_n3lo*(convSPa + convSPb)/static_cast<double>(_nf);

					// q
					for (uint j=1; j<=_nf; j++)
						HFT_N3LO1(arr[j], arr[j+6], j, k, SP);
					// qbar
					for (uint j=1+6; j<=_nf+6; j++)
						HFT_N3LO2(arr[j-6], arr[j], j, k, SP);

					HFT_N3LO3(arr_singlet[0], arr_singlet[1], k); // gluon
					HFT_N3LO4(arr_singlet[0], arr_singlet[1], k); // heavy flavor
				}
            } break;
		}
    }
    void DGLAPSolver::HFT_NNLO1(ArrayGrid& c, uint j, uint k)
    {
        double const as = _alpha_s.post(_nf+1);
        double const conv = getSplitFunc("A2ns").convolution(c, k);
        
		_C2[j][0][0][0][k] += std::pow(as/(4.0*PI), 2) * conv;
    }
    void DGLAPSolver::HFT_NNLO2(ArrayGrid& s1, ArrayGrid& s2, uint k)
    {
        double const as = _alpha_s.post(_nf+1);
        double const conv1 = getSplitFunc("A2gq").convolution(s1, k);
        double const conv2 = getSplitFunc("A2gg").convolution(s2, k);

		_S2[0][0][0][k] += std::pow(as/(4.0*PI), 2) * (conv1 + conv2);
    }
    void DGLAPSolver::HFT_NNLO3(ArrayGrid& s1, ArrayGrid& s2, uint k)
    {
        double const as = _alpha_s.post(_nf+1);
        double const conv1 = getSplitFunc("A2hq").convolution(s1, k);
        double const conv2 = getSplitFunc("A2hg").convolution(s2, k);
		double const fac = 0.5*std::pow(as/(4.0*PI), 2) * (conv1 + conv2);

		_C2[_nf+1][0][0][0][k] = fac;
        _C2[_nf+7][0][0][0][k] = fac;
    }

    // q
	void DGLAPSolver::HFT_N3LO1(ArrayGrid& q, ArrayGrid& qb, uint j, uint k, double SP)
	{
	    const double as = _alpha_s.post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

		double const conv1a = getSplitFunc("A2ns").convolution(q, k);
		double const conv1b = getSplitFunc("A3nsp").convolution(q, k);

		double const conv2a = getSplitFunc("A2ns").convolution(qb, k);
		double const conv2b = getSplitFunc("A3nsp").convolution(qb, k);

		double const conv3a = conv1a;
		double const conv3b = getSplitFunc("A3nsm").convolution(q, k);

		double const conv4a = conv2a;
		double const conv4b = getSplitFunc("A3nsm").convolution(qb, k);

		_D2[j][0][0][0][0][k] += 0.5*(
			((fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b)) +
			((fac_nnlo*conv3a + fac_n3lo*conv3b) - (fac_nnlo*conv4a + fac_n3lo*conv4b)) +
			SP);
	}

    // qbar
	void DGLAPSolver::HFT_N3LO2(ArrayGrid& q, ArrayGrid& qb, uint j, uint k, double SP)
	{
		const double as = _alpha_s.post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

		double const conv1a = getSplitFunc("A2ns").convolution(q, k);
		double const conv1b = getSplitFunc("A3nsp").convolution(q, k);

		double const conv2a = getSplitFunc("A2ns").convolution(qb, k);
		double const conv2b = getSplitFunc("A3nsp").convolution(qb, k);

		double const conv3a = conv1a;
		double const conv3b = getSplitFunc("A3nsm").convolution(q, k);

		double const conv4a = conv2a;
		double const conv4b = getSplitFunc("A3nsm").convolution(qb, k);

		_D2[j][0][0][0][0][k] += 0.5*(
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
	    
		const double conv1a = getSplitFunc("A2gq").convolution(qp, k);
		const double conv1b = getSplitFunc("A3gq").convolution(qp, k);
		const double conv2a = getSplitFunc("A2gg").convolution(g, k);
		const double conv2b = getSplitFunc("A3gg").convolution(g, k);
		
		_S2[0][0][0][k] += (fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b);
	}

	// heavy flavor
	void DGLAPSolver::HFT_N3LO4(ArrayGrid& g, ArrayGrid& qp, uint k)
	{
		const double as = _alpha_s.post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

	    const double conv1a = getSplitFunc("A2hq").convolution(qp, k);
		const double conv1b = getSplitFunc("A3hq").convolution(qp, k);
		const double conv2a = getSplitFunc("A2hg").convolution(g, k);
		const double conv2b = getSplitFunc("A3hg").convolution(g, k);
		
        const double res = 0.5*((fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b));
		_D2[(_nf+1)][0][0][0][0][k] = res;
        _D2[(_nf+1)+6][0][0][0][0][k] = res;
	}
} // namespace Candia2
