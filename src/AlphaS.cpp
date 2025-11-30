#include "Candia-v2/AlphaS.hpp"
#include "Candia-v2/Common.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <numbers>

namespace Candia2
{
	// helper function to ensure that NF is not out of bounds
	static void __assert_nf(const uint nf)
	{
		if (nf > 8)
		{
			std::cerr << "[ERROR] AlphaS::__assert_nf(): Encountered bad nf value: " << nf << '\n';
			exit(1);
		}
	}
	
	double AlphaS::Masses(const uint nf) const
	{
		__assert_nf(nf);
		return _masses[nf];
	}

	double AlphaS::Pre(const uint nf) const
	{
		__assert_nf(nf);
		return _pre[nf];
	}

	double AlphaS::Post(const uint nf) const
	{
		__assert_nf(nf);
		return _post[nf];
	}

	

	
	double AlphaS::CalcBeta0(const uint nf) const
	{
		return (11.0/3.0)*NC - (4.0/3.0)*TR*static_cast<double>(nf);
	}

	double AlphaS::CalcBeta1(const uint nf) const
	{
	    return (34.0/3.0)*NC*NC - (4.0*CF + (20.0/3.0)*NC)*TR*static_cast<double>(nf);
	}

	double AlphaS::CalcBeta2(const uint nf) const
	{
		double f = static_cast<double>(nf);
		return (2857.0/54.0)*NC*NC*NC + (2*CF*CF - (205.0/9.0)*CF*NC - (1415.0/27.0)*NC*NC)*TR*f +
			((44.0/9.0)*CF + (158.0/27.0)*NC)*TR*TR*f*f;
	}

	double AlphaS::CalcBeta3(const uint nf) const
	{
		const double F = static_cast<double>(nf);
		return (149753.0/6.0 + 3564.0*Zeta3)
			   - (1078361.0/162.0 + (6508.0/27.0)*Zeta3)*F
			   + (50065.0/162.0 + (6472.0/81.0)*Zeta3)*F*F
			   + (1093.0/729.0)*F*F*F;
	}

	double AlphaS::BetaFn(const double alpha) const
	{
		double res = _beta0;

		if (_order > 0)
			res += _beta1*alpha/(4.0*PI);
		if (_order > 1)
			res += _beta2*std::pow(alpha,2)/(16.0*PI_2);
		if (_order > 2)
			res += _beta3*std::pow(alpha,3)/(64.0*PI_3);
		
		res *= -std::pow(alpha,2)/(4.0*PI);

		return res;
	}


	double AlphaS::PreMatch(const double alpha, const uint nf)
	{
		double Nf = static_cast<double>(nf);
		if (_order == 0)
			return alpha;

		double a = alpha;
		double res = alpha;

		if (_order >= 1)
			res += -a*a*(1.0/6.0)*_L/PI;
		if (_order >= 2)
			res += a*a*a*((1.0/36.0)*_L*_L - (19.0/24.0)*_L - (7.0/24.0))/PI_2;
		if (_order >= 3)
		{
		    res += std::pow(alpha, 4)*(
				-(58933.0/124416.0) - (2.0/3.0)*Zeta2 - (2.0/9.0)*Zeta2*std::numbers::log2e - (80507.0/27648.0)*Zeta3
				- (8521.0/1728.0)*_L - (131.0/576.0)*_L*_L - (1.0/216.0)*_L*_L*_L
				+ Nf*((2479.0/31104.0) + (1.0/9.0)*Zeta2 + (409.0/1728.0)*_L)
			)/PI_3;
		}

		return res;
	}

	double AlphaS::PostMatch(const double alpha, const uint nf)
	{
		double Nf = static_cast<double>(nf);
		if (_order == 0)
			return alpha;

		double a = alpha;
		double res = alpha;

	    if (_order >= 1)
			res += a*a*(1.0/6.0)*_L/PI;
		if (_order >= 2)
			res += a*a*a*((7.0/24.0) + (19.0/24.0)*_L + (1.0/36.0)*_L*_L)/PI_2;
		if (_order >= 3)
		{
			res += std::pow(alpha, 4)*(
				(58933.0/124416.0) + (2.0/3.0)*Zeta2 + (2.0/9.0)*Zeta2*std::numbers::log2e + (80507.0/27648.0)*Zeta3
				+ (8521.0/1728.0)*_L + (131.0/576.0)*_L*_L + (1.0/216.0)*_L*_L*_L
				- Nf*((2479.0/31104.0) + (1.0/9.0)*Zeta2 + (409.0/1728.0)*_L)
			)/PI_3;
		}

		return res;
	}


	void AlphaS::CalculateThresholdValues(const double Qf)
	{
		uint nfi = 3,
			 nff = 6;
		uint nf1;
		const double mu0 = std::sqrt(2.0);

		// here we fix the final value in our masses/energy array
		// to correspond to the final energy that we are evolving to
		double temp = Qf;
		for (nff=6; temp<=_masses[nff]; nff--);

		uint idx;
		if (temp>_masses[6])
			idx = 7;
		else
			for (idx=nfi+1; temp>_masses[idx]; idx++);

		_masses[idx] = temp;
		for (uint i=idx+1; i<8; i++)  // TODO: why does this loop inclusive on 8?
			_masses[i] = 0.0;

		
		// set nf1 correctly
		for (nf1=nff; mu0<_masses[nf1]; nf1--);
		if (nf1<nfi)
			nf1++;

		std::cerr << std::fixed << std::setprecision(9);
		std::cerr << "[ALPHAS] initial: nf1 = " << nf1 << "\talpha0 = " << _alpha0 << '\n';

		Update(nf1);

		_post[nf1] = Evaluate(mu0, _masses[nf1], _alpha0);
		_pre[nf1]  = PreMatch(_post[nf1], nf1);

		_pre[nf1+1] =  Evaluate(mu0, _masses[nf1+1], _alpha0);
		_post[nf1+1] = PostMatch(_pre[nf1+1], nf1);

		uint nf;
		for (nf=nf1-1; nf>=nfi; nf--)
		{
			Update(nf);
			
			_post[nf] = Evaluate(_masses[nf+1], _masses[nf], _pre[nf+1]);
			_pre[nf]  = PreMatch(_post[nf], nf);
		}

		for (nf=nf1+1; nf<=nff+1; nf++)
		{
			Update(nf-1);
			
			_pre[nf]  = Evaluate(_masses[nf-1], _masses[nf], _post[nf-1]);
			_post[nf] = PostMatch(_pre[nf], nf);
		}

		std::cerr << "[ALPHAS] Computed alpha_s threshold values. They are:\n";

		std::cerr << std::fixed << std::setprecision(9);
		for (nf=nfi; nf<=nff+1; nf++)
			std::cerr << "[ALPHAS] " << nf << ' '
					  << std::setw(14) << _masses[nf] << ' '
					  << std::setw(14) << _pre[nf] << ' '
					  << std::setw(14) << _post[nf] << '\n';

	}


	double AlphaS::Evaluate(const double Qi, const double Qf, const double alpha0) const
	{
		// if the before/after energies are identical, there is nothing to evaluate
		if (Qi == Qf)
			return alpha0;

		if (Qf < Qi)
		{
			std::cerr << "[ERROR] AlphaS::Evaluate(): Final energy Qf=" << Qf << " is smaller than initial energy Qi=" << Qi
					  << ". Will return alpha0=" << alpha0 << std::endl;
			return alpha0;
		}

		// at LO we have the exact solution
		if (_order == 0) {
			return (2.0*PI*alpha0) / (2.0*PI + alpha0*_beta0*log(Qf/Qi));
		}

		// otherwise, 4th order runge-kutta
		const static uint steps = 200;
		double h = 2.0*std::log(Qf/Qi) / static_cast<double>(steps);
		double k1, k2, k3, k4;
		double a = alpha0;
		
		for (uint i=0; i<steps; i++) {
			k1 = h*BetaFn(a);
			k2 = h*BetaFn(a + k1/2.0);
			k3 = h*BetaFn(a + k2/2.0);
			k4 = h*BetaFn(a + k3);
			
			a += (k1/6.0) + (k2/3.0) + (k3/3.0) + (k4/6.0);
		}

		return a;
	}




	void AlphaS::Update(const uint nf)
	{
		_beta0 = CalcBeta0(nf);
		_beta1 = CalcBeta1(nf);
		_beta2 = CalcBeta2(nf);
		_beta3 = CalcBeta3(nf);
	}





	uint AlphaS::Nff(const uint nfi, const double Qf)
	{
	    double aux = Qf;
		uint nff, i;
		
		for (nff=6; aux<=_masses[nff]; nff--);

		if (aux>_masses[6])
			i=7;
		else
			for (i=nfi+1; aux>_masses[i]; i++);

		_masses[i] = aux;
		for (uint j=i+1; j<8; j++)
			_masses[j]=0.;

		return nff;
	}

}
