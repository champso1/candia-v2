#include "Candia-v2/SplittingFn.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/AlphaS.hpp"
#include "Candia-v2/SpecialFuncs.hpp"

namespace Candia2
{		
	uint SplittingFunction::_nf = 3;
	
	double P0ns::Regular(const double x) const
	{
		return CF*(-1.0-x);
	}
	double P0ns::Plus(const double x) const
	{
		UNUSED(x);
		return 2.0*CF;
	}
	double P0ns::Delta(const double x) const
	{
		UNUSED(x);
		return (3.0/2.0)*CF;
	}


	double P0qq::Regular(const double x) const
	{
		return CF*(-1.0-x);
	}
	double P0qq::Plus(const double x) const
	{
		UNUSED(x);
		return 2.0*CF;
	}
	double P0qq::Delta(const double x) const
	{
		UNUSED(x);
		return (3.0/2.0)*CF;
	}


	double P0qg::Regular(const double x) const
	{
		return 2.0*TR*_nf*(2.0*x*x - 2.0*x + 1.0);
	}
	double P0qg::Plus(const double x) const
	{
		UNUSED(x);
		return 0;
	}
	double P0qg::Delta(const double x) const
	{
		UNUSED(x);
		return 0;
	}


	double P0gq::Regular(const double x) const
	{
		return CF*(x - 2.0 + 2.0/x);
	}
	double P0gq::Plus(const double x) const
	{
		UNUSED(x);
		return 0;
	}
	double P0gq::Delta(const double x) const
	{
		UNUSED(x);
		return 0;
	}


	double P0gg::Regular(const double x) const
	{
		return 2.0*NC*(1.0/x - 2.0 + x - x*x);
	}
	double P0gg::Plus(const double x) const
	{
		UNUSED(x);
		return 2.0*NC;
	}
	double P0gg::Delta(const double x) const
	{
		UNUSED(x);
		return _beta0/2.0;
	}


	double P1nsp::Regular(const double x) const
	{
		double f1 = CF*(4.0*_nf*TR*(1.0 + x)*(-1.0 + 11.0*x) + NC*(89.0 + (-134.0 + 6.0*M_PI_2 - 223.0*x)*x) + 6.0*CF*(-27.0 + M_PI_2 + (27.0 + M_PI_2)*x*x));
		f1 /= 18.0*(1.0 + x);
		
		double f2 = Li2(-x)*(2.0*CF*(2.0*CF - NC)*(1.0 + x*x));
		f2 /= 1.0+x;
		
		double f3 = std::log(x)*(CF*(30.0*CF - 230.*NC + 4.0*_nf*TR + 12.0*CF*x + (-24.0*CF + NC + 4.0*_nf*TR)*x*x));
		f3 /= 6.0*(-1.0 + x);
		
		double f4 = std::pow(std::log(x), 2.0)*(CF*(2.0*NC*(1.0 + x*x) + CF*(-1.0 + x)*(3.0 + x*(2.0 + 3.0*x))));
		f4 /= 2.0*(-1.0 + x*x);

		double f5 = std::log(x)*std::log1p(-x)*(2.0*CF*CF*(1.0 + x*x));
		f5 /= -1.0+x;

		double f6 = std::log(x)*std::log1p(x)*(2.0*CF*(2.0*CF - NC)*(1.0 + x*x));
		f6 /= 1.0+x;
		return f1 + f2 + f3 - f4 + f5 + f6;
	}
	double P1nsp::Plus(const double x) const
	{
		UNUSED(x);
		return -(CF*(NC*(-67.0 + 3.0*M_PI_2) + 20.0*_nf*TR))/9.0;
	}
	double P1nsp::Delta(const double x) const
	{
		UNUSED(x);
		return (CF*(-4.0*_nf*(3.0 + 4.0*M_PI_2)*TR + NC*(51.0 + 44.0*M_PI_2 - 216.0*Zeta3) + 9.0*CF*(3.0 - 4.0*M_PI_2 + 48.0*Zeta3)))/72.0;
	}

	
	double P1nsm::Regular(const double x) const
	{
		double xp1 = 1.0+x;
		
		double f1 = (CF*(4.0*_nf*TR*xp1*(-1.0 + 11.0*x) - 6.0*CF*(3.0 + M_PI_2 + (-3.0 + M_PI_2)*x*x) + NC*(-(xp1*(-17.0 + 151.0*x)) + 6.0*M_PI_2*(1 + x + x*x))));
		f1 /= 18.0*xp1;

		double f2 = Li2(-x)*(-2.0*CF*(2.0*CF - NC)*(1.0 + x*x));
		f2 /= xp1;

		double f3 = std::log(x)*(CF*(6.0*CF*(1.0 + 2.0*x) - (11.0*NC - 4.0*_nf*TR)*(1.0 + x*x)));
		f3 /= 6.0*(-1.0 + x);

		double f4 = std::pow(std::log(x), 2.0)*(CF*(CF*std::pow(-1.0+x, 3.0) - 2.0*NC*x*(1.0 + x*x)));
		f4 /= 2.0*(-1.0 + x*x);

		double f5 = std::log(x)*std::log1p(-x)*(2.0*CF*CF*(1.0 + x*x));
		f5 /= -1.0+x;

		double f6 = std::log(x)*std::log1p(x)*(-2.0*CF*(2.0*CF - NC)*(1.0 + x*x));
		f6 /= xp1;
		
		return f1 + f2 + f3 + f4 + f5 + f6;
	}
	double P1nsm::Plus(const double x) const
	{
		UNUSED(x);
		return -(CF*(NC*(-67.0 + 3.0*M_PI_2) + 20.0*_nf*TR))/9.0;
	}
	double P1nsm::Delta(const double x) const
	{
		UNUSED(x);
		return (CF*(-4.0*_nf*(3.0 + 4.0*M_PI_2)*TR + NC*(51.0+44.0*M_PI_2 - 216.0*Zeta3) + 9.0*CF*(3.0 - 4.0*M_PI_2 + 48.0*Zeta3)))/72.0;
	}


	double P1qq::Regular(const double x) const
	{
		double xp1 = 1.0+x;
		
		double f1 = CF*(4.0*_nf*TR*(20.0 + x + 46.0*x*x + 9.0*std::pow(x, 3.0) - 56.0*std::pow(x, 4.0)) + x*(-6.0*CF*(3.0 + M_PI_2 + (-3.0 + M_PI_2)*x*x) + NC*(-(xp1*(-17.0 + 151.0*x)) + 6.0*M_PI_2*(xp1 + x*x))));
		f1 /= 18.0*x*xp1;

		double f2 = Li2(-x)*(-2.0*CF*(2.0*CF - NC)*(1.0 + x*x));
		f2 /= xp1;

		double f3 = std::log(x)*(CF*(6.0*CF*(1.0 + 2.0*x) - 11.0*NC*(1.0 + x*x) + 8.0*_nf*TR*(-1.0 + 2.0*x*(-3.0 + 2.0*x*xp1))));
		f3 /= 6.0*(-1.0 + x);

		double f4 = std::pow(std::log(x), 2.0)*(CF*(CF*std::pow(-1.0+x, 3.0) - 2.0*(2.0*_nf*TR*(-1.0 + x)*pow(1.0+x, 2.0) + NC*x*(1.0 + x*x))));
		f4 /= 2.0*(-1.0 + x*x);

		double f5 = std::log(x)*std::log1p(-x)*(2.0*CF*CF*(1.0 + x*x));
		f5 /= -1.0 + x;

		double f6 = std::log(x)*std::log1p(x)*(-2.0*CF*(2.0*CF - NC)*(1.0 + x*x));
		f6 /= xp1;
		
		return f1 + f2 + f3 + f4 + f5 + f6;
	}
	double P1qq::Plus(const double x) const
	{
		UNUSED(x);
		return -(CF*(NC*(-67.0 + 3.0*M_PI_2) + 20.0*_nf*TR))/9.0;
	}
	double P1qq::Delta(const double x) const
	{
		UNUSED(x);
		return (CF*(-4.0*_nf*(3.0 + 4.0*M_PI_2)*TR + NC*(51.0 + 44.0*M_PI_2 - 216.0*Zeta3) + 9.0*CF*(3.0 - 4.0*M_PI_2 + 48.0*Zeta3)))/72.0;
	}


	double P1qg::Regular(const double x) const
	{
		double f1 = _nf*TR*(3.0*CF*x*(42.0 - 87.0*x + 60.0*x*x + M_PI_2*(-2.0 - 4.0*(-1.0 + x)*x)) - 2.0*NC*(-20.0 + x*(18.0 + x*(-225.0 + 6.0*M_PI_2 + 218.0*x))));
		f1 /= 9.0*x;

		double f2 = Li2(-x)*4.0*NC*_nf*TR*(1.0 + 2.0*x*(1.0 + x));
		double f3 = std::log(x)*(_nf*TR*(6.0*NC + 8.0*NC*x*(6.0 + 11.0*x) + 3.0*CF*(3.0 - 4.0*x + 8.0*x*x)))/3.0;
		double f4 = std::log1p(-x)*8.0*(CF-NC)*_nf*TR*(-1.+x)*x;
		double f5 = std::pow(std::log(x), 2.0)*_nf*TR*(CF - 2.0*NC - 2.0*(CF + 2.0*NC)*x + 4.0*CF*x*x);
		double f6 = std::pow(std::log1p(-x), 2.0)*2.0*(CF-NC)*_nf*TR*(1.0 + 2.0*(-1.+x)*x);
		double f7 = std::log(x)*std::log1p(-x)*4.0*CF*_nf*TR*(1.0 + 2.0*(-1.0+x)*x);
		double f8 = std::log(x)*std::log1p(x)*4.0*NC*_nf*TR*(1.0 + 2.0*x*(1.0+x));
		
		return f1 - f2 + f3 - f4 + f5 + f6 - f7 - f8;
	}
	double P1qg::Plus(const double x) const
	{
		UNUSED(x);
		return 0;
	}
	double P1qg::Delta(const double x) const
	{
		UNUSED(x);
		return 0;
	}


	double P1gq::Regular(const double x) const
	{
		double f1 = CF*(-9.0*CF*x*(5.0 + 7.0*x) - 16.0*_nf*TR*(5.0 + x*(-5.0 + 4.0*x)) + 2.0*NC*(9.0 + x*(19.0 + 6.0*M_PI_2 + x*(37.0 + 44.0*x))));
		f1 /= 18.0*x;

		double f2 = Li2(-x)*(2.0*CF*NC*(2.0 + x*(2.0+x)))/x;
		double f3 = std::log(x)*(CF*(3.0*CF*(4.0 + 7.0*x) - 2.0*NC*(36.0 + x*(15.0 + 8.0*x))))/6.0;
		
		double f4 = std::log1p(-x)*(CF*(-4.0*_nf*TR*(2.0 + (-2.0+x)*x) - 3.0*CF*(6.0 + x*(-6.0 + 5.0*x)) + NC*(22.0 + x*(-22.0 + 17.0*x))));
		f4 /= 3.0*x;

		double f5 = std::pow(std::log(x), 2.0)*(CF*(CF*(-2.0+x) + 2.0*NC*(2.0+x)))/2.0;
		double f6 = std::pow(std::log1p(-x), 2.0)*((CF*(CF-NC)*(2.0 + (-2.0+x)*x))/x);
		double f7 = std::log(x)*std::log1p(-x)*(-2.0*CF*NC*(2.0 + (-2.+x)*x))/x;
		double f8 = std::log(x)*std::log1p(x)*(2.0*CF*NC*(2.0 + x*(2.0+x)))/x;
		
		return f1 + f2 + f3 + f4 + f5 - f6 + f7 + f8;
	}
	double P1gq::Plus(const double x) const
	{
		UNUSED(x);
		return 0;
	}
	double P1gq::Delta(const double x) const
	{
		UNUSED(x);
		return 0;
	}


	double P1gg::Regular(const double x) const
	{
		double xp1 = 1.0+x;

		double f1 = 24.0*CF*_nf*TR*(-1.0+x)*xp1*(-1.0 + x*(11.0 + 5.0*x)) + NC*(NC*x*(-(xp1*(25.0 + 109.0*x)) + 6.0*M_PI_2*(3.0 + 2.0*x*(2.0 + x + x*x))) + 4.0*_nf*TR*(-23.0 + x*(6.0 + x*(10.0 + x*(4.0 + 23.0*x)))));
		f1 /= 18.0*x*xp1;

		double f2 = Li2(-x)*(4.0*NC*NC*std::pow(xp1 + x*x, 2.0));
		f2 /= x*xp1;

		double f3 = std::log(x)*(-4.0*NC*_nf*TR*xp1 - 6.0*CF*_nf*TR*(3.0 + 5.0*x) + NC*NC*(-25.0 + 11.0*(1.0 - 4.0*x)*x))/3.0;

		double f4 = std::pow(std::log(x), 2.0)*(-2.0*(CF*_nf*TR*(-1.0+x)*pow(xp1, 2.0) + NC*NC*std::pow(-1.0 + (-1.0+x)*x, 2.0)));
		f4 /= -1.0 + x*x;

		double f5 = std::log(x)*std::log1p(-x)*(4.0*NC*NC*std::pow(1.0 + (-1.0+x)*x, 2.0));
		f5 /= (-1.0+x)*x;

		double f6 = std::log(x)*std::log1p(+x)*(4.0*NC*NC*std::pow(xp1 + x*x, 2.0));
		f6 /= x*xp1;
		
		return f1 + f2 + f3 + f4 + f5 + f6;
	}
	double P1gg::Plus(const double x) const
	{
		UNUSED(x);
		return -(NC*(NC*(-67.0 + 3.0*M_PI_2) + 20.0*_nf*TR))/9.0;
	}
	double P1gg::Delta(const double x) const
	{
		UNUSED(x);
		return -(CF*_nf*TR) + (NC*(-4.0*_nf*TR + NC*(8.0 + 9.0*Zeta3)))/3.0;
	}



	double P2nsp::Regular(const double x) const
	{
		const double dl = std::log(x);
		const double dl1 = std::log1p(-x);
		const double d81 = 1.0/81.0;

		const double nf = static_cast<double>(_nf);

		double res = 1641.1 - 3135.0*x + 243.6*std::pow(x, 2.0) - 522.1*std::pow(x, 3.0)
			+ 128.*d81*std::pow(dl, 4.0) + 2400.*d81*std::pow(dl, 3.0)
			+ 294.9*std::pow(dl, 2.0) + 1258.0*dl
			+ 714.1*dl1 + dl*dl1*(563.9 + 256.8*dl)
			+ nf * ( -197.0 + 381.1*x + 72.94*std::pow(x, 2.0) + 44.79*std::pow(x, 3.0)
					 - 192.0*d81*std::pow(dl, 3.0) - 2608.0*d81*std::pow(dl, 2.0) - 152.6*dl
					 - 5120.0*d81*dl1 - 56.66*dl*dl1 - 1.497*x*std::pow(dl, 3.0) )
			+ std::pow(nf, 2.0)*( 32.0*x*dl/(1.0-x) * (3.0*dl + 10.0) + 64.0
								  + (48.0*std::pow(dl, 2.0) + 352.0*dl + 384.0)*(1.0-x) ) * d81;

		return res/8.0;
	}
	double P2nsp::Plus(const double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		
		double res = (1174.898 - nf*183.187 - pow(nf, 2)*(64.0/81.0));
		
		return res/8.0;
	}
	double P2nsp::Delta(const double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		// const double dl1 = std::log1p(-x);

		double res = 1295.624 - 0.24 - nf*(173.938 - 0.011) + std::pow(nf, 2)*1.13067;
		
		return res/8.0;
	}


	double P2nsm::Regular(const double x) const
	{
		const double dl = std::log(x);
		const double dl1 = std::log1p(-x);
		const double d81 = 1.0/81.0;

		const double nf = static_cast<double>(_nf);

		double res = 1860.2 - 3505.0*x + 297.0*std::pow(x, 2.0) - 433.2*std::pow(x, 3.0)
			+ 116.0*d81*std::pow(dl, 4.0) + 2880.0*d81*std::pow(dl, 3.0) 
			+ 399.2*std::pow(dl, 2.0) + 1465.2*dl
			+ 714.1*dl1 + dl*dl1*(684.0 + 251.2*dl)
			+ nf*( -216.62 + 406.5*x + 77.89*std::pow(x, 2.0) + 34.76*std::pow(x, 3.0)
				   - 256.0*d81*std::pow(dl, 3.0) - 3216.0*d81*std::pow(dl, 2.0) - 172.69*dl 
				   - 5120.*d81*dl1 - 65.43*dl*dl1 - 1.136*x*std::pow(dl, 3.0) )
			+ std::pow(nf, 2.0)*( 32.0*x*dl/(1.0-x) * (3.0*dl + 10.0) + 64.0
								  + (48.0*std::pow(dl, 2.0) + 352.0*dl + 384.0)*(1.0-x) ) * d81;

		return res/8.0;
	}
	double P2nsm::Plus(const double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		
		double res = (1174.898 - nf*183.187 - pow(nf, 2)*(64.0/81.0));
		
		return res/8.0;
	}
	double P2nsm::Delta(const double x) const
	{
		UNUSED(x);

		const double nf = static_cast<double>(_nf);

		double res = (1295.624 - 0.154) - nf*(173.938 - 0.005) + std::pow(nf, 2.0)*(1.13067);

		return res/8.0;
	}


	double P2nsv::Regular(const double x) const
	{
		const double dl = std::log(x);
		const double x1 = 1.0-x;
		const double dl1 = std::log1p(-x);
		const double d27 = 1.0/27.0;
		const double d81 = 1.0/81.0;
		
		const double nf = static_cast<double>(_nf);

		double res1 = 1860.2 - 3505.0*x + 297.0*std::pow(x, 2.0) - 433.2*std::pow(x, 3.0)
			+ 116.0*d81*std::pow(dl, 4.0) + 2880.0*d81*std::pow(dl, 3.0) 
			+ 399.2*std::pow(dl, 2.0) + 1465.2*dl
			+ 714.1*dl1 + dl*dl1*(684.0 + 251.2*dl)
			+ nf*( -216.62 + 406.5*x + 77.89*std::pow(x, 2.0) + 34.76*std::pow(x, 3.0)
				   - 256.0*d81*std::pow(dl, 3.0) - 3216.0*d81*std::pow(dl, 2.0) - 172.69*dl 
				   - 5120.*d81*dl1 - 65.43*dl*dl1 - 1.136*x*std::pow(dl, 3.0) )
			+ std::pow(nf, 2.0)*( 32.0*x*dl/(1.0-x) * (3.0*dl + 10.0) + 64.0
								  + (48.0*std::pow(dl, 2.0) + 352.0*dl + 384.0)*(1.0-x) ) * d81;

		double res2 = x1* ( 151.49 + 44.51 * x - 43.12 * std::pow(x, 2) + 4.820 * std::pow(x, 3) )
			+ 40.*d27 * std::pow(dl, 4) - 80.*d27 * std::pow(dl, 3) + 6.892 * std::pow(dl, 2) 
			+ 178.04 * dl + dl*dl1 * ( - 173.1 + 46.18 * dl )
			+ x1*dl1 * ( - 163.9 / x - 7.208 * x );
		res2 *= nf;

		return (res1 + res2)/8.0;
	}
	double P2nsv::Plus(const double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		
		double res = (1174.898 - nf*183.187 - pow(nf, 2)*(64.0/81.0));
		
		return res/8.0;
	}
	double P2nsv::Delta(const double x) const
	{
		UNUSED(x);

		const double nf = static_cast<double>(_nf);

		double res = (1295.624 - 0.154) - nf*(173.938 - 0.005) + std::pow(nf, 2.0)*(1.13067);

		return res/8.0;
	}


	double P2qq::Regular(const double x) const
	{
		const double dl = std::log(x);
		const double dl1 = std::log1p(-x);
		const double d81 = 1.0/81.0;

		const double nf = static_cast<double>(_nf);

		double res1 = 1641.1 - 3135.0*x + 243.6*std::pow(x, 2.0) - 522.1*std::pow(x, 3.0)
			+ 128.*d81*std::pow(dl, 4.0) + 2400.*d81*std::pow(dl, 3.0)
			+ 294.9*std::pow(dl, 2.0) + 1258.0*dl
			+ 714.1*dl1 + dl*dl1*(563.9 + 256.8*dl)
			+ nf * ( -197.0 + 381.1*x + 72.94*std::pow(x, 2.0) + 44.79*std::pow(x, 3.0)
					 - 192.0*d81*std::pow(dl, 3.0) - 2608.0*d81*std::pow(dl, 2.0) - 152.6*dl
					 - 5120.0*d81*dl1 - 56.66*dl*dl1 - 1.497*x*std::pow(dl, 3.0) )
			+ std::pow(nf, 2.0)*( 32.0*x*dl/(1.0-x) * (3.0*dl + 10.0) + 64.0
								  + (48.0*std::pow(dl, 2.0) + 352.0*dl + 384.0)*(1.0-x) ) * d81;

		double res2a = 3584.0/(27.0*x) * dl - 506.0/x + 160.0/27.0 * std::pow(dl, 4.0)
			- 400.0/9.0 * std::pow(dl, 3.0) + 131.4*std::pow(dl, 2.0) - 661.6*dl
			- 5.926*std::pow(dl1, 3.0) - 9.751*std::pow(dl1, 2.0) - 72.11*dl1
			+ 177.4 + 392.9*x - 101.4*std::pow(x, 2.0) - 57.04*dl*dl1;
		double res2b = 256.0/(81.0*x) + 32.0/27.0 * std::pow(dl, 3.0) + 17.89*std::pow(dl, 2.0)
			+ 61.75*dl + 1.778*std::pow(dl1, 2.0) + 5.944*dl1 + 100.1
			- 125.2*x + 49.26*std::pow(x, 2.0) - 12.59*std::pow(x, 3.0)
			- 1.889*dl*dl1;

		double res2 = (1.0-x)*nf*(res2a + nf*res2b);

		return (res1 + res2)/8.0;
	}
	double P2qq::Plus(const double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		
		double res = (1174.898 - nf*183.187 - pow(nf, 2)*(64.0/81.0));
		
		return res/8.0;
	}
	double P2qq::Delta(const double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		// const double dl1 = std::log1p(-x);

		double res = 1295.624 - 0.24 - nf*(173.938 - 0.011) + std::pow(nf, 2)*1.13067;
		
		return res/8.0;
	}


	double P2qg::Regular(const double x) const
	{
		const double nf = static_cast<double>(_nf);

		const double dl = std::log(x);
		const double dl1 = std::log1p(-x);

		
		double res1 = - 896.0/(3.0*x) * dl - 1268.3/x + 536.0/27.0 * std::pow(dl, 4.0) 
			- 44.0/3.0 * std::pow(dl, 3.0) + 881.5*std::pow(dl, 2.0) + 424.9*dl 
			+ 100.0/27.0 * std::pow(dl1, 4.0) - 70.0/9.0 * std::pow(dl1, 3.0) 
			- 120.5*std::pow(dl1, 2.0) + 104.42*dl1
			+ 2522.0 - 3316.0*x + 2126.0*std::pow(x, 2.0)
			+ dl*dl1*(1823.0 - 25.22*dl) - 252.5*x*std::pow(dl, 3.0);

		
		double res2 = 1112.0/(243.0*x) - 16.0/9.0 * std::pow(dl, 4.0) 
			- 376.0/27.0 * std::pow(dl, 3.0) - 90.8*std::pow(dl, 2.0) - 254.0*dl   
			+ 20.0/27.0 * std::pow(dl1, 3.0) + 200.0/27.0 * std::pow(dl1, 2.0) - 5.496*dl1
			- 252.0  + 158.0*x + 145.4*std::pow(x, 2.0) - 139.28*std::pow(x, 3.0)
			- dl*dl1*(53.09 + 80.616*dl) - 98.07*x*std::pow(dl, 2.0)
			+ 11.70*x*std::pow(dl, 3.0);

		double res = nf*(res1 + nf*res2);
		return res/8.0;
	}
	double P2qg::Plus(const double x) const
	{
		UNUSED(x);
		return 0;
	}
	double P2qg::Delta(const double x) const
	{
		UNUSED(x);
		return 0;
	}


	double P2gq::Regular(const double x) const
	{
		const double nf = static_cast<double>(_nf);

		const double dl = std::log(x);
		const double dl1 = std::log1p(-x);

		
		double res1 = 1189.3 * dl/x + 6163.1/x - 4288.0/81.0 * std::pow(dl, 4.0)
			+ 1568.0/9.0 * std::pow(dl, 3.0) - 1794.0*std::pow(dl, 2.0) + 4033.0*dl
			+ 400.0/81.0 * std::pow(dl1, 4.0) + 2200.0/27.0 * std::pow(dl1, 3.0)
			+ 606.3*std::pow(dl1, 2.0) + 2193.0*dl1 
			- 4307.0 + 489.3*x + 1452.0*std::pow(x, 2.0) + 146.0*std::pow(x, 3.0)
			- 447.3*std::pow(dl, 2.0)*dl1 - 972.9*x*std::pow(dl, 2.0);
		double res2 = 71.082 * dl/x - 46.41/x + 128.0/27.0 * std::pow(dl, 4.0)
			+ 704.0/81.0 * std::pow(dl, 3.0) + 20.39*std::pow(dl, 2.0) + 174.8*dl
			- 400.0/81.0 * std::pow(dl1, 3.0) - 68.069*std::pow(dl1, 2.0) - 296.7*dl1
			- 183.8 + 33.35*x - 277.9*std::pow(x, 2.0) + 108.6*x*std::pow(dl, 2.0)
			- 49.68*dl*dl1;
		double res3 = 64.0*(-1.0/x + 1.0 + 2.0*x)
			+ 320.0*dl1*(1.0/x - 1.0 + 0.8*x)
			+ 96.0*std::pow(dl1, 2.0)*(1.0/x - 1.0 + 0.5*x);
		res3 /= 27.0;

		double res = res1 + nf*(res2 + nf*res3);
		return res/8.0;
	}
	double P2gq::Plus(const double x) const
	{
		UNUSED(x);
		return 0;
	}
	double P2gq::Delta(const double x) const
	{
		UNUSED(x);
		return 0;
	}


	double P2gg::Regular(const double x) const
	{
		const double nf = static_cast<double>(_nf);

		const double dl = std::log(x);
		const double dl1 = std::log1p(-x);

		
		double res1 = 2675.8 * dl/x + 14214.0/x - 144.0*std::pow(dl, 4.0) + 72.0*std::pow(dl, 3.0)
			- 7471.0*std::pow(dl, 2.0) + 274.4*dl + 3589.0*dl1 - 20852.0
			+ 3968.0*x - 3363.0*std::pow(x, 2.0) + 4848.0*std::pow(x, 3.0) 
			+ dl*dl1*(7305.0 + 8757.0*dl);
		double res2 = 157.27 * dl/x + 182.96/x + 512.0/27.0 * std::pow(dl, 4.0)
			+ 832.0/9.0 * std::pow(dl, 3.0) + 491.3*std::pow(dl, 2.0) + 1541.0*dl
			- 320.0*dl1 - 350.2 + 755.7*x - 713.8*std::pow(x, 2.0) 
			+ 559.3*std::pow(x, 3.0) + dl*dl1*(26.15 - 808.7*dl);
		double res3 = -680.0/(243.0*x) - 32.0/27.0 * std::pow(dl, 3.0) + 9.680*std::pow(dl, 2.0)
			- 3.422*dl - 13.878 + 153.4*x - 187.7*std::pow(x, 2.0) 
			+ 52.75*std::pow(x, 3.0) - dl*dl1*(115.6 - 85.25*x + 63.23*dl);

		double res = res1 + nf*(res2 + nf*res3);
		return res/8.0;
	}
	double P2gg::Plus(const double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		
		double res = 2643.521 - nf*412.172 - std::pow(nf, 2.0)*(16.0/9.0);

		return res/8.0;
	}
	double P2gg::Delta(const double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);

		//const double dl1 = std::log1p(-x);

		double res = (4425.448 + 0.446) - nf*(528.720 + 0.003) + std::pow(nf, 2.0)*(6.4630);
		
		return res/8.0;
	}




	
	double P3nsp::Regular(const double x) const
	{
		
		const double x2   = x * x;
		const double x3   = x2 * x;
		const double omx  = 1 - x;
		const double dm   = 1 / omx;
		const double dl   = log(x);
		const double dl2  = dl * dl;
		const double dl3  = dl2 * dl;
		const double dl4  = dl3 * dl;
		const double dl5  = dl4 * dl;
		const double dl6  = dl5 * dl;
		const double dlm  = log(omx);
		const double dlm2 = dlm * dlm;
		const double dlm3 = dlm2 * dlm;

		// Leading large-n_c, nf^0 and nf^1, parametrized
		const double p3nsa0  =
			2.5e+4 * ( omx * ( 3.5254 + 8.6935 * x - 1.5051 * x2 + 1.8300 * x3 )
					   + 11.883 * x * dl - 0.09066 * x * dl2 + 11.410 * omx * dlm + 13.376  * dl * dlm )
			+ 5.167133e+4 * dl + 1.712095e+4 * dl2 + 2.863226e+3 * dl3 + 2.978255e+2 * dl4
			+ 1.6e+1 * dl5 + 5.e-1 * dl6 - 2.973385e+4 + 1.906980e+4 * dlm;
		const double p3nsa1  =
			2.5e+4 * ( omx * ( - 0.74077 + 1.4860 * x - 0.23631 * x2 + 0.31584 * x3 )
					   + 2.5251 * omx * dlm + 2.5203 * dl * dlm + 2.2242 * x * dl
					   - 0.02460 * x * dl2 + 0.00310 * x * dl3 )
			- 9.239374e+3 * dl - 2.917312e+3 * dl2 - 4.305308e+2 *dl3 - 3.6e+1 * dl4
			- 4. / 3. * dl5 + 8.115605e+3 - 3.079761e+3 * dlm;

		// Nonleading large-n_c, nf^0 and nf^1: two approximations
		const double p3npa01 =
			3948.16 * omx - 2464.61 * ( 2 * x - x2 ) * omx - 1839.44 * dl2 - 402.156 * dl3
			- 1777.27 * dlm2 * omx - 204.183 * dlm3 * omx + 507.152 - 5.587553e1 * dl4 - 2.831276e0 * dl5
			- 1.488340e-1 * dl6 - 2.601749e3 - 2.118867e3 * dlm;
		const double p3npa02 =
			( 8698.39 - 10490.47 * x ) * x * omx + 1389.73 * dl + 189.576 * dl2
			- 173.936 * dlm2 * omx + 223.078 * dlm3 * omx + 505.209 - 5.587553e1 * dl4 - 2.831276e0 * dl5
			- 1.488340e-1 * dl6 - 2.601749e3 - 2.118867e3 * dlm;

		const double p3npa11 =
			( - 1116.34 + 1071.24 * x ) * x * omx - 59.3041 * dl2 - 8.4620 * dl3
			- 143.813 * dlm * omx - 18.8803 * dlm3 * omx - 7.33927 + 4.658436e0*dl4 + 2.798354e-1 * dl5
			+ 3.121643e2 + 3.379310e2 * dlm;
		const double p3npa12 =
			( - 690.151 - 656.386 * x2 ) * omx + 133.702 * dl2 + 34.0569 * dl3
			- 745.573 * dlm * omx + 8.61438 * dlm3 * omx - 7.53662 + 4.658437e0 * dl4 + 2.798354e-1 * dl5
			+ 3.121643e2 + 3.379310e2 * dlm;

		// nf^2 (parametrized) and nf^3 (exact)
		const double p3nspa2 =
			2.5e+2 * ( omx * ( 3.0008 + 0.8619 * x - 0.12411 * x2 + 0.31595* x3 )
- 0.37529 * x * dl - 0.21684 * x * dl2 - 0.02295 * x * dl3
					   + 0.03394 * omx * dlm + 0.40431 * dl * dlm )
			+ 3.930056e+2 * dl + 1.125705e+2 * dl2 + 1.652675e+1 * dl3
			+ 7.901235e-1 * dl4 - 3.760092e+2 + 2.668861e+1 * dlm;
		const double p3nsa3  =
			- 2.426296 - 8.460488e-1 * x + ( 5.267490e-1 * dm - 3.687243 + 3.160494 * x ) * dl
			- ( 1.316872 * ( dm + 1e-1) - 1.448560 * x ) * dl2
			- ( 2.633745e-1 * dm - 1.31687e-1 * ( 1 + x ) ) * dl3;

		// Assembly
		const double p3nspai = p3nsa0 + _nf * p3nsa1 + _nf * _nf * p3nspa2 + _nf * _nf * _nf * p3nsa3;
		if (_imod == 1)
			return p3nspai + p3npa01 + _nf * p3npa11;
		else if (_imod == 2)
			return p3nspai + p3npa02 + _nf * p3npa12;
		else
			return p3nspai + 0.5 * ( ( p3npa01 + p3npa02 ) + _nf * ( p3npa11 + p3npa12 ) );
	}
	double P3nsp::Plus(const double x) const
	{
		const double d1 = 1 / ( 1 - x );

		const double a4qi =
			2.120902e+4
			- 5.179372e+3 * _nf
			+ 1.955772e+2 * _nf * _nf
			+ 3.272344e+0 * _nf * _nf * _nf;
		const double a4ap1 = - 507.152 + 7.33927 * _nf;
		const double a4ap2 = - 505.209 + 7.53662 * _nf;

		if (_imod == 1)
			return ( a4qi + a4ap1 ) * d1;
		else if (_imod == 2)
			return ( a4qi + a4ap2 ) * d1;
		else
			return ( a4qi + 0.5 * ( a4ap1 + a4ap2 ) ) * d1;
	}
	double P3nsp::Delta(const double x) const
	{
		const double dl1 = log(1 - x);

		const double a4qi  =
			2.120902e+4
			- 5.179372e+3 * _nf
			+ 1.955772e+2 * _nf * _nf
			+ 3.272344e+0 * _nf * _nf * _nf;
		const double a4ap1 = - 507.152 + 7.33927 * _nf;
		const double a4ap2 = - 505.209 + 7.53662 * _nf;

		const double b4qi =
			2.579609e+4 + 0.08
			- ( 5.818637e+3 + 0.97 )   * _nf
			+ ( 1.938554e+2 + 0.0037 ) * _nf * _nf
			+   3.014982e+0 * _nf * _nf * _nf;
		const double b4ap1 = - 2405.03 + 267.965 * _nf;
		const double b4ap2 = - 2394.47 + 269.028 * _nf;

		if (_imod == 1)
			return ( a4qi + a4ap1 ) * dl1 + b4qi + b4ap1;
		else if (_imod == 2)
			return ( a4qi + a4ap1 ) * dl1 + b4qi + b4ap2;
		else
			return ( a4qi + 0.5 * ( a4ap1 + a4ap2 ) ) * dl1 + b4qi + 0.5 * ( b4ap1 + b4ap2 );
	}


	double P3nsm::Regular(const double x) const
	{
		const double x2   = x * x;
		const double x3   = x2 * x;
		const double omx  = 1 - x;
		const double dm   = 1 / omx;
		const double dl   = log(x);
		const double dl2  = dl * dl;
		const double dl3  = dl2 * dl;
		const double dl4  = dl3 * dl;
		const double dl5  = dl4 * dl;
		const double dl6  = dl5 * dl;
		const double dlm  = log(omx);
		const double dlm2 = dlm * dlm;
		const double dlm3 = dlm2 * dlm;

		// Leading large-n_c, nf^0 and nf^1, parametrized
		const double p3nsa0  =
			2.5e+4 * ( omx * ( 3.5254 + 8.6935 * x - 1.5051 * x2 + 1.8300 * x3 )
					   + 11.883 * x * dl - 0.09066 * x * dl2 + 11.410 * omx * dlm + 13.376  * dl * dlm )
			+ 5.167133e+4 * dl + 1.712095e+4 * dl2 + 2.863226e+3 * dl3 + 2.978255e+2 * dl4
			+ 1.6e+1 * dl5 + 5.e-1 * dl6 - 2.973385e+4 + 1.906980e+4 * dlm;
		const double p3nsa1  =
			2.5e+4 * ( omx * ( - 0.74077 + 1.4860 * x - 0.23631 * x2 + 0.31584 * x3 )
					   + 2.5251 * omx * dlm + 2.5203 * dl * dlm + 2.2242 * x * dl
					   - 0.02460 * x * dl2 + 0.00310 * x * dl3 )
			- 9.239374e+3 * dl - 2.917312e+3 * dl2 - 4.305308e+2 *dl3 - 3.6e+1 * dl4
			- 4. / 3. * dl5 + 8.115605e+3 - 3.079761e+3 * dlm;

		// Nonleading large-n_c, nf^0 and nf^1: two approximations
		const double p3nma01 =
			( 5992.88 * ( 1 + 2 * x ) + 31321.44 * x2 ) * omx + 511.228 - 1618.07 * dl + 2.25480 * dl3
			+ 31897.82 * dlm * omx + 4653.76 * dlm2 * omx + 4.964335e-1 * ( dl6 + 6 * dl5 )
			- 2.601749e+3 - 2.118867e+3 * dlm;
		const double p3nma02 =
			( 4043.59 - 15386.6 * x ) * x * omx + 502.481 + 1532.96  * dl2 + 31.6023 * dl3
			- 3997.39  * dlm * omx + 511.567 * dlm3 * omx + 4.964335e-1 * ( dl6 + 18 * dl5 )
			- 2.601749e+3 - 2.118867e+3 * dlm;

		const double p3nma11 =
			( 114.457 * ( 1 + 2 * x ) + 2570.73 * x2 ) * omx - 7.08645 - 127.012 * dl2 + 2.69618 * dl4
			+ 1856.63 * dlm * omx + 440.17 * dlm2 * omx + 3.121643e+2 + 3.379310e+2 * dlm;
		const double p3nma12 =
			( - 335.995 * ( 2 + x ) - 1605.91 * x2 ) * omx - 7.82077 - 9.76627 * dl2 + 0.14218 * dl5
			- 1360.04 * dlm * omx + 38.7337 * dlm3 * omx + 3.121643e+2 + 3.379310e+2 * dlm;

		// nf^2 (parametrized) and nf^3 (exact)
		const double p3nsma2 =
			2.5e+2 * ( omx * ( 3.2206 + 1.7507 * x + 0.13281 * x2 + 0.45969 * x3 )
					   + 1.5641 * x * dl - 0.37902 * x * dl2 - 0.03248 * x *dl3
					   + 2.7511 * omx * dlm + 3.2709  * dl * dlm )
			+ 4.378810e+2 * dl + 1.282948e+2 * dl2 + 1.959945e+1 * dl3
			+ 9.876543e-1 * dl4 - 3.760092e+2 + 2.668861e+1 * dlm;
		const double p3nsa3  =
			- 2.426296 - 8.460488e-1 * x + ( 5.267490e-1 * dm - 3.687243 + 3.160494 * x ) * dl
			- ( 1.316872 * ( dm + 1e-1) - 1.448560 * x ) * dl2
			- ( 2.633744e-1 * dm - 1.31687e-1 * ( 1 + x ) ) * dl3;

		// Assemblx
		const double p3nsmai = p3nsa0 + _nf * p3nsa1 + _nf * _nf * p3nsma2 + _nf * _nf * _nf * p3nsa3;
		if (_imod == 1)
			return p3nsmai + p3nma01 + _nf * p3nma11;
		else if (_imod == 2)
			return p3nsmai + p3nma02 + _nf * p3nma12;
		else
			return p3nsmai + 0.5 * ( ( p3nma01 + p3nma02 ) + _nf * ( p3nma11 + p3nma12 ) );
	}
	double P3nsm::Plus(const double x) const
	{
		const double d1 = 1 / ( 1 - x );

		const double a4qi  =
			2.120902e+4
			- 5.179372e+3 * _nf
			+ 1.955772e+2 * _nf * _nf
			+ 3.272344e+0 * _nf * _nf * _nf;
		const double a4ap1 = - 511.228 + 7.08645 * _nf;
		const double a4ap2 = - 502.481 + 7.82077 * _nf;

		if (_imod == 1)
			return ( a4qi + a4ap1 ) * d1;
		else if (_imod == 2)
			return ( a4qi + a4ap2 ) * d1;
		else
			return ( a4qi + 0.5 * ( a4ap1 + a4ap2 ) ) * d1;
	}
	double P3nsm::Delta(const double x) const
	{
		const double dl1 = log(1 - x);

		const double a4qi  =
			2.120902e+4
			- 5.179372e+3 * _nf
			+ 1.955772e+2 * _nf * _nf
			+ 3.272344e+0 * _nf * _nf * _nf;
		const double a4ap1 = - 511.228 + 7.08645 * _nf;
		const double a4ap2 = - 502.481 + 7.82077 * _nf;

		const double b4qi =
			2.579609e+4 + 0.08
			- ( 5.818637e+3 + 0.97 )   * _nf
			+ ( 1.938554e+2 + 0.0037 ) * _nf * _nf
			+   3.014982e+0 * _nf * _nf * _nf;
		const double b4ap1 = - 2426.05  + 266.674 * _nf - 0.05 * _nf;
		const double b4ap2 = - 2380.255 + 270.518 * _nf - 0.05 * _nf;

		if (_imod == 1)
			return ( a4qi + a4ap1 ) * dl1 + b4qi + b4ap1;
		else if (_imod == 2)
			return ( a4qi + a4ap1 ) * dl1 + b4qi + b4ap2;
		else
			return ( a4qi + 0.5 * ( a4ap1 + a4ap2 ) ) * dl1 + b4qi + 0.5 * ( b4ap1 + b4ap2 );
	}


	double P3nss::Regular(const double x) const
	{
		const double x2   = x * x;
		const double omx  = 1 - x;
		const double dl   = log(x);
		const double dl2  = dl * dl;
		const double dl3  = dl2 * dl;
		const double dl4  = dl3 * dl;
		const double dl5  = dl4 * dl;
		const double dl6  = dl5 * dl;
		const double dlm  = log(omx);
		const double dlm2 = dlm * dlm;
		const double dlm3 = dlm2 * dlm;

		// nf^1: two approximations
		const double p3nsa11 =
			omx * x * ( 4989.2 - 1607.73 * x ) + 3687.6 * dl + 3296.6 * dl2 + 1271.11* dl3
			+ 533.44 * dl4 + 97.27 *  dl5 + 4 * dl6 + 60.40 * omx * dlm2 + 4.685 * omx * dlm3;
		const double p3nsa12 =
			1030.79 * omx * x + 1266.77 * omx * ( 2 - x2 ) + 2987.83 * dl + 273.05 * dl2 - 923.48 * dl3
			- 236.76 * dl4 - 33.886 * dl5 - 4 * dl6 - 254.63 * omx * dlm - 0.28953 * omx * dlm3;

		// nf^2 (parametrized)
		const double p3nssa2 =
			2.5e+2 * ( omx * ( - 4.7656 + 1.6908 * x + 0.1703 * x2 )
					   - 0.41652 * x *dl + 0.90777 * x * dl2 + 0.12478 * x * dl3
					   + 0.17155 * omx * dlm + 0.17191  * dl * dlm )
			- 6.473971e+2 * dl - 6.641219e+1 * dl2 - 5.353347 * dl3 - 5.925926 * dl4
			- 3.950617e-1 * dl5 + 1.970002e+1 * omx * dlm - 3.435474 * omx * dlm2;

		if (_imod == 1)
			return _nf * p3nsa11 + _nf * _nf * p3nssa2;
		else if (_imod == 2)
			return _nf * p3nsa12 + _nf * _nf * p3nssa2;
		else
			return 0.5 *_nf * ( p3nsa11 + p3nsa12 ) + _nf * _nf * p3nssa2;
	}
	double P3nss::Plus(const double x) const
	{
		UNUSED(x);
		return 0.0;
	}
	double P3nss::Delta(const double x) const
	{
		UNUSED(x);
		return 0.0;
	}

} // namespace Candia2
