#include "Candia-v2/Candia.hpp"

#include <iostream>

namespace Candia2
{

	void DGLAPSolver::HeavyFlavorTreatmentAlt()
	{
		if (_order < 2)
		{
			std::cerr << "[DGLAP] HeavyFlavorTreatmentAlt(): Order is not NNLO or higher. Skipping heavy flavor treatment...\n";
			return;
		}
		
		std::cerr << "[DGLAP] ALT: Treating heavy flavors: "
				  << _nf+1
				  << "th quark mass threshold (mass "
				  << _alpha_s.Masses(_nf+1) << ")\n";

		OpMatElem::update(-_log_mur2_muf2, _nf);

		// Copy of pre-threshold distributions
		// the nf+1 dists are defined in terms of the nf dists,
		// so we need this copy since we would otherwise be overwriting
		// a given dist as we do the calculations
		// we just store them in the s=1 array, and modify the s=0 array
		// (which are the initial conditions for the next set of iterations
		// at the next nf
		std::cerr << "[DGLAP] Creating copy of pre-threshold distributions... ";

		for (uint j=0; j<=1; ++j)
			std::ranges::copy(_S[0][j][0], _S[0][j][1].begin());
		
		if (_order == 2)
		{
			for (uint i=1; i<=_nf; i++)
			{
				for (uint j=i; j<=i+6; j+=6)
					std::ranges::copy(_C[j][0][0][0], _C[j][1][0][0].begin());
			}
			
		}
		else if (_order == 3)
		{
			for (uint i=1; i<=_nf; i++)
			{
				for (uint j=i; j<=i+6; j+=6)
					std::ranges::copy(_D[j][0][0][0][0], _D[j][1][0][0][0].begin());
			}
		}
		std::cerr << "Done.\n";

		double as = _alpha_s.Post(_nf+1);
		std::cerr << "[DGLAP] Value of alpha_s post threshold: " << as << '\n';

		switch(_order)
		{
			case 2:
			{
				for (uint k=0; k<_grid.Size()-1;k++)
				{
					// q
					for (uint j=1; j<=_nf; j++)
						HFT_NNLO1_ORIG(j, k);
					// qbar
					for (uint j=1+6; j<=_nf+6; j++)
						HFT_NNLO1_ORIG(j, k);

					HFT_NNLO2_ORIG(k); // gluon
					HFT_NNLO3_ORIG(k); // heavy flavor
				}
			} break;
			case 3:
			{
				for (uint k=0; k<_grid.Size()-1;k++)
				{
					const double fac_n3lo = as*as*as/(64.0*PI_3);
					const double convSPa = _grid.Convolution(_S[0][1][1], _A3psqq, k);
					const double convSPb = _grid.Convolution(_S[0][0][1], _A3sqg, k);
					const double SP = fac_n3lo*(convSPa + convSPb)/static_cast<double>(_nf);


					for (uint j=1; j<=_nf; j++)
						HFT_N3LO1_ORIG(j, k);
					// qbar
					for (uint j=1+6; j<=_nf+6; j++)
						HFT_N3LO1_ORIG(j, k);

					HFT_N3LO2_ORIG(k); // gluon
					HFT_N3LO3_ORIG(k); // heavy flavor

					
					continue;
		
				    // q
					for (uint j=1; j<=_nf; j++)
					{
						HFT_N3LO1(j, k, SP);
					}
					// qbar
					for (uint j=1+6; j<=_nf+6; j++)
						HFT_N3LO2(j, k, SP);

					HFT_N3LO3(k); // gluon
					HFT_N3LO4(k); // heavy flavor
				}
			} break;
		}
	}
	
	
	void DGLAPSolver::HFT_NNLO1_ORIG(uint j, uint k)
	{
		const double as = _alpha_s.Post(_nf+1);
		_C[j][0][0][0][k] += std::pow(as/(4.0*PI), 2) * _grid.Convolution(_C[j][1][0][0], _A2ns, k);
	}
	void DGLAPSolver::HFT_NNLO2_ORIG(uint k)
	{
		const double as = _alpha_s.Post(_nf+1);
		_S[0][0][0][k] += std::pow(as/(4.0*PI), 2) * (_grid.Convolution(_S[0][1][1], _A2gq, k)
			+ _grid.Convolution(_S[0][0][1], _A2gg, k));
	}
	void DGLAPSolver::HFT_NNLO3_ORIG(uint k)
	{
		const double as = _alpha_s.Post(_nf+1);
		double fac = 0.5*std::pow(as/(4.0*PI), 2) * (_grid.Convolution(_S[0][1][1], _A2hq, k)
			+ _grid.Convolution(_S[0][0][1], _A2hg, k));
		_C[_nf+1][0][0][0][k] = _C[_nf+7][0][0][0][k] = fac;
	}


	void DGLAPSolver::HFT_N3LO1_ORIG(uint j, uint k)
	{
		const double as = _alpha_s.Post(_nf+1);
		_D[j][0][0][0][0][k] += std::pow(as/(4.0*PI), 2) * _grid.Convolution(_D[j][1][0][0][0], _A2ns, k);
	}
	void DGLAPSolver::HFT_N3LO2_ORIG(uint k)
	{
		const double as = _alpha_s.Post(_nf+1);
		_S[0][0][0][k] += std::pow(as/(4.0*PI), 2) * (_grid.Convolution(_S[0][1][1], _A2gq, k)
			+ _grid.Convolution(_S[0][0][1], _A2gg, k));
	}
	void DGLAPSolver::HFT_N3LO3_ORIG(uint k)
	{
		const double as = _alpha_s.Post(_nf+1);
		double fac = 0.5*std::pow(as/(4.0*PI), 2) * (_grid.Convolution(_S[0][1][1], _A2hq, k)
			+ _grid.Convolution(_S[0][0][1], _A2hg, k));
		_D[_nf+1][0][0][0][0][k] = _D[_nf+7][0][0][0][0][k] = fac;
	}


	// q
	void DGLAPSolver::HFT_N3LO1(uint j, uint k, double SP)
	{
	    const double as = _alpha_s.Post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

		const double conv1a = _grid.Convolution(_D[j][1][0][0][0], _A2ns, k);
		const double conv1b = _grid.Convolution(_D[j][1][0][0][0], _A3nsp, k);
		const double conv1c = conv1a;
		const double conv1d = _grid.Convolution(_D[j][1][0][0][0], _A3nsm, k);

		const double conv2a = _grid.Convolution(_D[j+6][1][0][0][0], _A2ns, k);
		const double conv2b = _grid.Convolution(_D[j+6][1][0][0][0], _A3nsp, k);
		const double conv2c = conv2a;
		const double conv2d = _grid.Convolution(_D[j+6][1][0][0][0], _A3nsm, k);
		
		_D[j][0][0][0][0][k] += 0.5*(
			((fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv1c + fac_n3lo*conv1d))
			+ ((fac_nnlo*conv2a + fac_n3lo*conv2b) - (fac_nnlo*conv2c + fac_n3lo*conv2d))
			+ SP);
	}

	// qbar
	void DGLAPSolver::HFT_N3LO2(uint j, uint k, double SP)
	{
		const double as = _alpha_s.Post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

		const double conv1a = _grid.Convolution(_D[j][1][0][0][0], _A2ns, k);
		const double conv1b = _grid.Convolution(_D[j][1][0][0][0], _A3nsp, k);
		const double conv1c = conv1a;
		const double conv1d = _grid.Convolution(_D[j][1][0][0][0], _A3nsm, k);

		const double conv2a = _grid.Convolution(_D[j+6][1][0][0][0], _A2ns, k);
		const double conv2b = _grid.Convolution(_D[j+6][1][0][0][0], _A3nsp, k);
		const double conv2c = conv2a;
		const double conv2d = _grid.Convolution(_D[j+6][1][0][0][0], _A3nsm, k);

		_D[j][0][0][0][0][k] += 0.5*(
			((fac_nnlo*conv1a + fac_n3lo*conv1b) - (fac_nnlo*conv1c + fac_n3lo*conv1d))
			+ ((fac_nnlo*conv2a + fac_n3lo*conv2b) + (fac_nnlo*conv2c + fac_n3lo*conv2d))
			+ SP);
	}

	// gluon (index 0 in S array)
	void DGLAPSolver::HFT_N3LO3(uint k)
	{
		const double as = _alpha_s.Post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);
	    
		const double conv1a = _grid.Convolution(_S[0][0][1], _A2gg, k);
		const double conv1b = _grid.Convolution(_S[0][0][1], _A3gg, k);
		const double conv2a = _grid.Convolution(_S[0][1][1], _A2gq, k);
		const double conv2b = _grid.Convolution(_S[0][1][1], _A3gq, k);
		
		_S[0][0][0][k] += (fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b);
	}

	// heavy flavor
	void DGLAPSolver::HFT_N3LO4(uint k)
	{
		const double as = _alpha_s.Post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

	    const double conv1a = _grid.Convolution(_S[0][1][1], _A2hq, k);
		const double conv1b = _grid.Convolution(_S[0][1][1], _A3hq, k);
		const double conv2a = _grid.Convolution(_S[0][0][1], _A2hg, k);
		const double conv2b = _grid.Convolution(_S[0][0][1], _A3hg, k);
		
		_D[_nf+1][0][0][0][0][k] = _D[_nf+7][0][0][0][0][k] = 0.5*((fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b));
	}
	
} // namespace Candia2
