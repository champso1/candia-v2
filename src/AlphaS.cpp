#include "Candia-v2/AlphaS.hpp"
#include "Candia-v2/Common.hpp"

#include <cstdlib>
#include <print>
#include <cmath>
#include <numbers>

namespace Candia2
{
	void AlphaS::assertNf() const
	{
		if (_nf > 8)
		{
			std::println("[AlphaS: ERROR] assertNf(): Found nf value of {}, expected < 8", _nf);
			exit(EXIT_FAILURE);
		}
	}

	void AlphaS::assertScheme() const
	{
		if (_scheme == UNSET)
		{
			std::println("[AlphaS: ERROR] assertScheme(): must set a scheme before accessing alpha_s or mass values.");
			exit(1);
		}
	}


	void AlphaS::setVFNS(std::array<double, 8> masses, uint nfi)
	{
		if (_scheme == FIXED)
		{
			std::println("[ALPHAS: WARNING] setVFNS(): scheme previously set to FFNS.");
		}
		else if (_scheme == VARIABLE)
		{
			std::println("[ALPHAS: WARNING] setVFNS(): scheme already set to VFNS. Overwriting previous masses...");
		}
	    _masses = masses;
		_scheme = VARIABLE;
		
		_nfi = nfi;
		double aux = _Qf;
		uint i{};
		
		for (_nff=6; aux<=_masses[_nff]; _nff--);

		if (aux>_masses[6])
			i=7;
		else
			for (i=nfi+1; aux>_masses[i]; i++);

		_masses[i] = aux;
		for (uint j=i+1; j<8; j++)
			_masses[j]=0.;

		calculateThresholdValues();
	}

	void AlphaS::setFFNS(uint nf)
	{
		if (_scheme == FIXED)
		{
			std::println("[ALPHAS: WARNING] setFFNS(): scheme already set to FFNS. Overwriting previous value of nf...");
		}
		else if (_scheme == VARIABLE)
		{
			std::println("[ALPHAS: WARNING] setFFNS(): scheme previously set to VFNS.");
		}
		_nf = nf;
		_nfi = nf;
		_nff = nf;
		_scheme = FIXED;

		// ensure the array is cleared
		_masses = ([]() -> std::array<double, 8>
		{
			std::array<double, 8> arr{};
			for (uint i=0; i<8; ++i)
				arr[i] = 0.0;
			return arr;
		})();

		_masses[_nfi] = _Q0;
		_masses[_nfi+1] = _Qf;

		// what we setup is the "post-threshold match" at the initial flavor
		// which is just the initial value of alpha_s, i.e. alpha0
		// and the "pre-threshold match" at the next flavor
		// which defines the final value of alpha_s,
		// this we just calculate.
		_post[_nfi] = _alpha0;
		_pre[_nfi+1] = evaluate(_Q0, _Qf, _alpha0);
	}

	
	double AlphaS::masses(uint nf) const
	{
		assertNf();
		assertScheme();
		return _masses[nf];
	}

	double AlphaS::pre(uint nf) const
	{
		assertNf();
		assertScheme();
		return _pre[nf];
	}

	double AlphaS::post(uint nf) const
	{
		assertNf();
		assertScheme();
		return _post[nf];
	}

	

	
	double AlphaS::calcBeta0(uint nf) const
	{
		return (11.0/3.0)*NC - (4.0/3.0)*TR*static_cast<double>(nf);
	}

	double AlphaS::calcBeta1(uint nf) const
	{
	    return (34.0/3.0)*NC*NC - (4.0*CF + (20.0/3.0)*NC)*TR*static_cast<double>(nf);
	}

	double AlphaS::calcBeta2(uint nf) const
	{
		double f = static_cast<double>(nf);
		return (2857.0/54.0)*NC*NC*NC + (2*CF*CF - (205.0/9.0)*CF*NC - (1415.0/27.0)*NC*NC)*TR*f +
			((44.0/9.0)*CF + (158.0/27.0)*NC)*TR*TR*f*f;
	}

	double AlphaS::calcBeta3(uint nf) const
	{
		const double F = static_cast<double>(nf);
		return (149753.0/6.0 + 3564.0*Zeta3)
			   - (1078361.0/162.0 + (6508.0/27.0)*Zeta3)*F
			   + (50065.0/162.0 + (6472.0/81.0)*Zeta3)*F*F
			   + (1093.0/729.0)*F*F*F;
	}

	double AlphaS::betaFn(double alpha) const
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


	double AlphaS::preMatch(double alpha, uint nf)
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

	double AlphaS::postMatch(double alpha, uint nf)
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


	void AlphaS::calculateThresholdValues()
	{
		// here we fix the final value in our masses/energy array
		// to correspond to the final energy that we are evolving to

		uint nf1{};
		// set nf1 correctly
		for (nf1=_nff; _Q0<_masses[nf1]; nf1--);
		if (nf1<_nfi)
			nf1++;

		std::println("[ALPHAS] initial: nf1 = {}\talpha0 = {}", nf1, _alpha0);

		update(nf1);

		_post[nf1] = evaluate(_Q0, _masses[nf1], _alpha0);
		_pre[nf1]  = preMatch(_post[nf1], nf1);

		_pre[nf1+1] =  evaluate(_Q0, _masses[nf1+1], _alpha0);
		_post[nf1+1] = postMatch(_pre[nf1+1], nf1);

		uint nf;
		for (nf=nf1-1; nf>=_nfi; nf--)
		{
			update(nf);
			
			_post[nf] = evaluate(_masses[nf+1], _masses[nf], _pre[nf+1]);
			_pre[nf]  = preMatch(_post[nf], nf);
		}

		for (nf=nf1+1; nf<=_nff+1; nf++)
		{
			update(nf-1);
			
			_pre[nf]  = evaluate(_masses[nf-1], _masses[nf], _post[nf-1]);
			_post[nf] = postMatch(_pre[nf], nf);
		}

		std::println("[ALPHAS] Computed alpha_s threshold values for VFNS. They are:");

		for (nf=_nfi; nf<=_nff+1; nf++)
			std::println("[ALPHAS] {} {:14.9} {:14.9} {:14.9}", nf, _masses[nf], _pre[nf], _post[nf]);

	}


	double AlphaS::evaluate(double Qi, double Qf, double alpha0) const
	{
		// if the before/after energies are identical, there is nothing to evaluate
		if (Qi == Qf)
			return alpha0;

		if (Qf < Qi)
		{
			std::println("[ALPHAS: ERROR] evaluate(): Final energy Qf={} is smaller than initial energy Qi={}.", Qf, Qi);
			exit(EXIT_FAILURE);
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
			k1 = h*betaFn(a);
			k2 = h*betaFn(a + k1/2.0);
			k3 = h*betaFn(a + k2/2.0);
			k4 = h*betaFn(a + k3);
			
			a += (k1/6.0) + (k2/3.0) + (k3/3.0) + (k4/6.0);
		}

		return a;
	}


	void AlphaS::update(uint nf)
	{
		_nf = nf;
		_beta0 = calcBeta0(nf);
		_beta1 = calcBeta1(nf);
		_beta2 = calcBeta2(nf);
		_beta3 = calcBeta3(nf);
	}

	


} // namespace Candia2
