#include "Candia-v2/AlphaS.hpp"

#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>

namespace Candia2
{	
	static void __assert_nf(const uint nf)
	{
		if (nf >= 8)
			throw("[AlphaS] Pre(): nf=" + std::to_string(nf) + " is out of bounds.");
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
		const double f = static_cast<double>(nf);
		
		return (149745.0/6.0 - 3564.0*Zeta3) - (1078361.0/162.0 + (6508.0/27.0)*Zeta3)*f
			+ (50065.0/162.0 + (6472.0/81.0)*Zeta3)*f*f + (1093.0/729.0)*f*f*f;
	}

	double AlphaS::BetaFn(const double alpha) const
	{
		double res = _beta0;

		if (_order > 0)
			res += _beta1*alpha/(4.0*M_PI);
		if (_order > 1)
			res += _beta2*std::pow(alpha,2)/(16.0*M_PI_2);
		if (_order > 2)
			res += _beta3*std::pow(alpha,3)/(64.0*M_PI_3);
		
		res *= -std::pow(alpha,2)/(4.0*M_PI);

		return res;
	}


	double AlphaS::PreMatch(const double alpha)
	{
		if (_order == 0)
			return alpha;
		
		double L = -std::log(1.0);
		double res = alpha;

		res += std::pow(alpha,2)*L / (6.0*M_PI);

		if (_order == 2)
			res += std::pow(alpha,3)*((L*L/36.0) - (19.0/24.0)*L - (7.0/24.0)) / (M_PI_2);

		return res;
	}

	double AlphaS::PostMatch(const double alpha)
	{
		if (_order == 0)
			return alpha;
		
		double L = -std::log(1.0);
		double res = alpha;

		res += std::pow(alpha,2)*L / (6.0*M_PI);

		if (_order == 2)
			res += std::pow(alpha,3)*(14.0 + 38.0*L + (4.0/3.0)*L*L) / (48.0*M_PI_2);

		return res;
	}


	void AlphaS::CalculateThresholdValues(const double Qf)
	{
		uint nfi = 3,
			 nff = 6;
		uint nf1;
		double mu0 = std::sqrt(2.0);

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
		for (int i=idx+1; i<=8; i++)
			_masses[i] = 0.0;

		
		std::cerr << "[ALPHA_S] Energy array: ";
		for (const double m : _masses)
		{
			std::cerr << m << ", ";
		}
		std::cerr << '\n';

		// set nf1 correctly
		for (nf1=nff; mu0<_masses[nf1]; nf1--);
		if (nf1<nfi)
			nf1++;

		std::cerr << "[ALPHA_S] nf1 = " << nf1 << "\talpha0 = " << _alpha0 << '\n';

		Update(nf1);

		_post[nf1] = Evaluate(mu0, _masses[nf1], _alpha0);
		_pre[nf1]  = PreMatch(_post[nf1]);

		_pre[nf1+1] =  Evaluate(mu0, _masses[nf1+1], _alpha0);
		_post[nf1+1] = PostMatch(_pre[nf1+1]);
		
		uint nf;
		for (nf=nf1-1; nf>=nfi; nf--)
		{
			Update(nf);
			_post[nf] = Evaluate(_masses[nf+1], _masses[nf], _pre[nf+1]);
			_pre[nf]  = PreMatch(_post[nf]);
		}

		
		for (nf=nf1+1; nf<=nff+1; nf++)
		{
			Update(nf-1);
			_pre[nf]  = Evaluate(_masses[nf-1], _masses[nf], _post[nf-1]);
			_post[nf] = PostMatch(_pre[nf]);
		}

		std::cerr << "[ALPHA_S] Computed alpha_s threshold values. They are:\n";
		
		for (nf=nfi; nf<=nff; nf++)
			std::cerr << "[ALPHA_S] " << nf << '\t' << _masses[nf] << '\t' << _pre[nf] << '\t' << _post[nf] << '\n';

	}


	double AlphaS::Evaluate(const double Qi, const double Qf, const double alpha0) const
	{
		if (Qi == Qf)
			return alpha0;

		if (_order == 0) {
			return (2.0*M_PI*alpha0) / (2.0F*M_PI + alpha0*_beta0*log(Qf/Qi));
		}

		const static uint steps = 200;
		double h = 2.0*std::log(Qf/Qi) / static_cast<double>(steps);
		double k1, k2, k3, k4;
		double a = alpha0;
		
		for (uint i=1; i<steps; i++) {
			k1 = h*BetaFn(a);
			k2 = h*BetaFn(a + k1/2.0);
			k3 = h*BetaFn(a + k2/2.0);
			k4 = h*BetaFn(a + k3);
			
			a += (k1/6.0) + (k2/3.0) + (k3/3.0) + (k4/6.0);
		}

		if (a > 1.0)
		{
			std::cerr << "[ALPHAS] Evaluate(): Error has occurred: Qi=" << Qi << ", Qf=" << Qf << ", alpha0=" << alpha0 << '\n';
			exit(1);
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

}
