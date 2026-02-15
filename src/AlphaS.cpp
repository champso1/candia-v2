#include "Candia-v2/AlphaS.hpp"

#include <cmath>

namespace Candia2
{
    void AlphaS::assertNf() const
	{
		if (_nf > 8)
			log(LOG_ERROR, "AlphaS::assertNf()", "Found nf value of {}, expected < 8", _nf);
	}

	void AlphaS::assertScheme() const
	{
		if (_scheme == UNSET)
			log(LOG_ERROR, "AlphaS::assertScheme()", "Must set a scheme before accessing alpha_s or mass values.");
	}

	void AlphaS::setVFNS(std::array<double, 8> const& masses, uint nfi)
	{
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

		log(LOG_DEBUG, "AlphaS", "Calculated mass array: {}", vec_to_str(_masses));
		calculateThresholdValues();
	}

	void AlphaS::setFFNS(uint nf)
	{
		_nf = nf;
		_nfi = nf;
		_nff = nf;
		_scheme = FIXED;

		// ensure the array is cleared
	    for (uint i=0; i<8; ++i)
			_masses[i] = 0.0;

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

		auto L = options.use_broken_log_value ? 1.0 : _L;
		auto L2 = L*L;
		auto L3 = L2*L;
		log(LOG_DEBUG, "AlphaS::preMatch()", "Using L={}", L);

		if (_order >= 1)
			res += -a*a*(1.0/6.0)*L/PI;
		if (_order >= 2)
			res += a*a*a*(-(7.0/24.0) - (19.0/24.0)*L + (1.0/36.0)*L2)/PI_2;
		if (_order >= 3) {
		    res += std::pow(alpha, 4)*(
				-(58933.0/124416.0) - (2.0/3.0)*Zeta2 - (2.0/9.0)*Zeta2*std::numbers::log2e - (80507.0/27648.0)*Zeta3
				- (8521.0/1728.0)*L - (131.0/576.0)*L2 - (1.0/216.0)*L3
				+ Nf*((2479.0/31104.0) + (1.0/9.0)*Zeta2 + (409.0/1728.0)*L))/PI_3;
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

		auto L = options.use_broken_log_value ? 1.0 : _L;
		auto L2 = L*L;
		auto L3 = L2*L;
		log(LOG_DEBUG, "AlphaS::postMatch()", "Using L={}", L);

	    if (_order >= 1)
			res += a*a*(1.0/6.0)*L/PI;
		if (_order >= 2)
			res += a*a*a*((7.0/24.0) + (19.0/24.0)*L + (1.0/36.0)*L2)/PI_2;
		if (_order >= 3) {
			res += std::pow(alpha, 4)*(
				(58933.0/124416.0) + (2.0/3.0)*Zeta2 + (2.0/9.0)*Zeta2*std::numbers::log2e + (80507.0/27648.0)*Zeta3
				+ (8521.0/1728.0)*L + (131.0/576.0)*L2 + (1.0/216.0)*L3
				- Nf*((2479.0/31104.0) + (1.0/9.0)*Zeta2 + (409.0/1728.0)*L))/PI_3;
		}

		return res;
	}


	void AlphaS::calculateThresholdValues()
	{
	    double mur_muf = std::sqrt(_mur2_muf2);
		log(LOG_DEBUG, "AlphaS::calculateThresholdValues()", "Using mur/muf={}", mur_muf);
		uint nf1 = 0;

		if (options.use_broken_log_value)
			log(LOG_WARNING, "AlphaS::calculateThresholdValues()", "Temporarily using broken value of the logarithm term (1.0).");
		
		for (nf1=_nff; _Q0<mur_muf*_masses[nf1]; nf1--);
		if (nf1<_nfi)
			nf1++;

		update(nf1);

		_post[nf1] = evaluate(_Q0, mur_muf*_masses[nf1], _alpha0);
		_pre[nf1]  = preMatch(_post[nf1], nf1); // might not need this value?

		_pre[nf1+1] =  evaluate(_Q0, mur_muf*_masses[nf1+1], _alpha0);
		_post[nf1+1] = postMatch(_pre[nf1+1], nf1);

		uint nf;
		for (nf=nf1-1; nf>=_nfi; nf--) {
			update(nf);
			
			_post[nf] = evaluate(mur_muf*_masses[nf+1], mur_muf*_masses[nf], _pre[nf+1]);
			_pre[nf]  = preMatch(_post[nf], nf);
		}

		for (nf=nf1+1; nf<=_nff+1; nf++) {
			update(nf-1);
			
			_pre[nf]  = evaluate(mur_muf*_masses[nf-1], mur_muf*_masses[nf], _post[nf-1]);
			_post[nf] = postMatch(_pre[nf], nf);
		}
		log(LOG_DEBUG, "AlphaS", "Computed alpha_s threshold values for VFNS. They are:");

		for (nf=_nfi; nf<=_nff+1; nf++)
			log(LOG_DEBUG, "AlphaS", "{} {:14.9} {:14.9} {:14.9}", nf, _masses[nf], _pre[nf], _post[nf]);

	}


	double AlphaS::evaluate(double Qi, double Qf, double alpha0) const
	{
		// if the before/after energies are identical, there is nothing to evaluate
		if (Qi == Qf)
			return alpha0;

		if (Qf < Qi)
			log(LOG_ERROR, "AlphaS::evaluate()", "Final energy Qf={} is smaller than initial energy Qi={}.", Qf, Qi);

		// at LO we have the exact solution
		if (_order == 0) {
			return (2.0*PI*alpha0) / (2.0*PI + alpha0*_beta0*std::log(Qf/Qi));
		}

		// otherwise, 4th order runge-kutta
		const static uint steps = 200;
		double h = 2.0*std::log(Qf/Qi) / static_cast<double>(steps);
		double k1{}, k2{}, k3{}, k4{};
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
