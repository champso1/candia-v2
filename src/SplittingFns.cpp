#include "Candia-v2/SplittingFn.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/SpecialFuncs.hpp"

#include <limits>
#include <cmath>

namespace Candia2
{		
	uint SplittingFunction::_nf = 4; // good default
	double SplittingFunction::_beta0 = 0.0; // will be set later
	double SplittingFunction::_log_muf2_mur2 = 0.0; // will be set later
	uint SplittingFunction::_imod = 3; //!< this is the avg
	
	double P0ns::calcRegular(double x) const
	{
		return CF*(-1.0-x);
	}
	double P0ns::calcPlus(double x) const
	{
		UNUSED(x);
		return 2.0*CF;
	}
	double P0ns::calcDelta(double x) const
	{
		UNUSED(x);
		return (3.0/2.0)*CF;
	}


	double P0qq::calcRegular(double x) const
	{
		return CF*(-1.0-x);
	}
	double P0qq::calcPlus(double x) const
	{
		UNUSED(x);
		return 2.0*CF;
	}
	double P0qq::calcDelta(double x) const
	{
		UNUSED(x);
		return (3.0/2.0)*CF;
	}


	double P0qg::calcRegular(double x) const
	{
		return 2.0*TR*_nf*(2.0*x*x - 2.0*x + 1.0);
	}


	double P0gq::calcRegular(double x) const
	{
		return CF*(x - 2.0 + 2.0/x);
	}


	double P0gg::calcRegular(double x) const
	{
		return 2.0*NC*(1.0/x - 2.0 + x - x*x);
	}
	double P0gg::calcPlus(double x) const
	{
		UNUSED(x);
		return 2.0*NC;
	}
	double P0gg::calcDelta(double x) const
	{
		UNUSED(x);
		return _beta0/2.0;
	}


	double P1nsp::calcRegular(double x) const
	{
		double Nf = static_cast<double>(_nf);
		return (CF*(4.*Nf*TR*(1.+x)*(-1.+11.*x)-6.*CF*(3.+PI_2+(-3.+PI_2)*x*x)+
				 NC*(-((1.+x)*(-17.+151.*x))+6.*PI_2*(1+x+x*x))))/(18.*(1+x))
				+Li2(-x)*(-2.*CF*(2.*CF-NC)*(1.+x*x))/(1.+x)
			    +std::log(x)*(CF*(6.*CF*(1.+2.*x)-(11.*NC-4.*Nf*TR)*(1.+x*x)))/(6.*(-1.+x))
			    +std::pow(std::log(x),2.)*(CF*(CF*std::pow(-1.+x,3.)-2.*NC*x*(1.+x*x)))/(2.*(-1.+x*x))
			    +std::log(x)*std::log(1.-x)*(2.*CF*CF*(1.+x*x))/(-1.+x)
			    +std::log(x)*std::log(1.+x)*(-2.*CF*(2.*CF-NC)*(1.+x*x))/(1.+x);
	}
	double P1nsp::calcPlus(double x) const
	{
		UNUSED(x);
		return -(CF*(NC*(-67.0 + 3.0*PI_2) + 20.0*_nf*TR))/9.0;
	}
	double P1nsp::calcDelta(double x) const
	{
		UNUSED(x);
		return (CF*(-4.0*_nf*(3.0 + 4.0*PI_2)*TR + NC*(51.0 + 44.0*PI_2 - 216.0*Zeta3) + 9.0*CF*(3.0 - 4.0*PI_2 + 48.0*Zeta3)))/72.0;
	}

	
	double P1nsm::calcRegular(double x) const
	{
		double Nf = static_cast<double>(_nf);
		return (CF*(4.*Nf*TR*(1.+x)*(-1.+11.*x)+NC*(89.+(-134.+6.*PI_2-223.*x)*x)+
				 6.*CF*(-27.+PI_2+(27.+ PI_2)*x*x)))/(18.*(1.+x))
				+Li2(-x)*(2.*CF*(2.*CF-NC)*(1.+x*x))/(1.+x)
			    +std::log(x)*(CF*(30.*CF-23.*NC+4.*Nf*TR+12.*CF*x+(-24.*CF+NC+4.*Nf*TR)*x*x))/(6.*(-1.+x))
			    -std::pow(std::log(x),2.)*(CF*(2.*NC*(1.+x*x)+CF*(-1.+x)*(3.+x*(2.+3.*x))))/(2.*(-1.+x*x))
			    +std::log(x)*std::log(1.-x)*(2.*CF*CF*(1.+x*x))/(-1.+x)
			    +std::log(x)*std::log(1.+x)*(2.*CF*(2.*CF-NC)*(1.+x*x))/(1.+x);
	}
	double P1nsm::calcPlus(double x) const
	{
		UNUSED(x);
		return -(CF*(NC*(-67.0 + 3.0*PI_2) + 20.0*_nf*TR))/9.0;
	}
	double P1nsm::calcDelta(double x) const
	{
		UNUSED(x);
		return (CF*(-4.0*_nf*(3.0 + 4.0*PI_2)*TR + NC*(51.0+44.0*PI_2 - 216.0*Zeta3) + 9.0*CF*(3.0 - 4.0*PI_2 + 48.0*Zeta3)))/72.0;
	}


	double P1qq::calcRegular(double x) const
	{
		double Nf = static_cast<double>(_nf);
		return (CF*(4.*Nf*TR*(20.+x+46.*x*x+9.*std::pow(x,3.)-56.*std::pow(x,4.))+
				     x*(-6.*CF*(3.+PI_2+(-3.+PI_2)*x*x)+NC*(-((1.+x)*(-17.+151.*x))+6.*PI_2*(1.+x+x*x)))))/(18.*x*(1.+x))
				+Li2(-x)*(-2.*CF*(2.*CF-NC)*(1.+x*x))/(1.+x)
			    +std::log(x)*(CF*(6.*CF*(1.+2.*x)-11.*NC*(1.+x*x)+8.*Nf*TR*(-1.+2.*x*(-3.+2.*x*(1.+x)))))/(6.*(-1.+x))
			    +std::pow(std::log(x),2.)*(CF*(CF*std::pow(-1.+x,3.)-2.*(2.*Nf*TR*(-1.+x)*std::pow(1.+x,2.)+NC*x*(1.+x*x))))/(2.*(-1.+x*x))
			    +std::log(x)*std::log(1.-x)*(2.*CF*CF*(1.+x*x))/(-1.+x)
			    +std::log(x)*std::log(1.+x)*(-2.*CF*(2.*CF-NC)*(1.+x*x))/(1.+x);
	}
	double P1qq::calcPlus(double x) const
	{
		UNUSED(x);
		return -(CF*(NC*(-67.0 + 3.0*PI_2) + 20.0*_nf*TR))/9.0;
	}
	double P1qq::calcDelta(double x) const
	{
		UNUSED(x);
		return (CF*(-4.0*_nf*(3.0 + 4.0*PI_2)*TR + NC*(51.0 + 44.0*PI_2 - 216.0*Zeta3) + 9.0*CF*(3.0 - 4.0*PI_2 + 48.0*Zeta3)))/72.0;
	}


	double P1qg::calcRegular(double x) const
	{
		double Nf = static_cast<double>(_nf);
		return (Nf*TR*(3.*CF*x*(42.-87.*x+60.*x*x+PI_2*(-2.-4.*(-1.+x)*x))-
				        2.*NC*(-20.+x*(18.+x*(-225.+6.*PI_2+218.*x)))))/(9.*x)
				-Li2(-x)*4.*NC*Nf*TR*(1.+2.*x*(1.+x))
			    +std::log(x)*(Nf*TR*(6.*NC+8.*NC*x*(6.+11.*x)+3.*CF*(3.-4.*x+8.*x*x)))/3.
			    -std::log(1.-x)*8.*(CF-NC)*Nf*TR*(-1.+x)*x
			    +std::pow(std::log(x),2.)*Nf*TR*(CF-2.*NC-2.*(CF+2.*NC)*x+4.*CF*x*x)
			    +std::pow(std::log(1.-x),2.)*2.*(CF-NC)*Nf*TR*(1.+2.*(-1.+x)*x)
			    -std::log(x)*std::log(1.-x)*4.*CF*Nf*TR*(1.+2.*(-1.+x)*x)
			    -std::log(x)*std::log(1.+x)*4.*NC*Nf*TR*(1.+2.*x*(1.+x));
	}


	double P1gq::calcRegular(double x) const
	{
		double Nf = static_cast<double>(_nf);
		return (CF*(-9.*CF*x*(5.+7.*x)-16.*Nf*TR*(5.+x*(-5.+4.*x))+2.*NC*(9.+x*(19.+6.*PI_2+x*(37.+44.*x)))))/(18.*x)
				+Li2(-x)*(2.*CF*NC*(2.+x*(2.+x)))/x
			    +std::log(x)*(CF*(3.*CF*(4.+7.*x)-2.*NC*(36.+x*(15.+8.*x))))/6.
			    +std::log(1.-x)*(CF*(-4.*Nf*TR*(2.+(-2.+x)*x)-3.*CF*(6.+x*(-6.+5.*x))+NC*(22.+x*(-22.+17.*x))))/(3.*x)
			    +std::pow(std::log(x),2.)*(CF*(CF*(-2.+x)+2.*NC*(2.+x)))/2.
			    -std::pow(std::log(1.-x),2.)*((CF*(CF-NC)*(2.+(-2.+x)*x))/x)
			    +std::log(x)*std::log(1.-x)*(-2.*CF*NC*(2.+(-2.+x)*x))/x
			    +std::log(x)*std::log(1.+x)*(2.*CF*NC*(2.+x*(2.+x)))/x;
	}



	double P1gg::calcRegular(double x) const
	{
		double Nf = static_cast<double>(_nf);
		return (24.*CF*Nf*TR*(-1.+x)*(1.+x)*(-1.+x*(11.+5.*x))+
				 NC*(NC*x*(-((1.+x)*(25.+109.*x))+6.*PI_2*(3.+2.*x*(2.+x+x*x)))+
				     4.*Nf*TR*(-23.+x*(6.+x*(10.+x*(4.+23.*x))))))/(18.*x*(1.+x))
			    +Li2(-x)*(4.*NC*NC*std::pow(1.+x+x*x,2.))/(x*(1.+x))
			    +std::log(x)*(-4.*NC*Nf*TR*(1.+x)-6.*CF*Nf*TR*(3.+5.*x)+NC*NC*(-25.+11.*(1.-4.*x)*x))/3.
			    +std::pow(std::log(x),2.)*(-2.*(CF*Nf*TR*(-1.+x)*std::pow(1+x,2.)+NC*NC*std::pow(-1.+(-1.+x)*x,2.)))/(-1.+x*x)
			    +std::log(x)*std::log(1.-x)*(4.*NC*NC*std::pow(1.+(-1.+x)*x,2.))/((-1.+x)*x)
			    +std::log(x)*std::log(1.+x)*(4.*NC*NC*std::pow(1.+x+x*x,2.))/(x*(1.+x));
	}
	double P1gg::calcPlus(double x) const
	{
		UNUSED(x);
		return -(NC*(NC*(-67.0 + 3.0*PI_2) + 20.0*_nf*TR))/9.0;
	}
	double P1gg::calcDelta(double x) const
	{
		UNUSED(x);
		return -(CF*_nf*TR) + (NC*(-4.0*_nf*TR + NC*(8.0 + 9.0*Zeta3)))/3.0;
	}



	double P2nsp::calcRegular(double x) const
	{
		const double dl = std::log(x);
		const double dlm = std::log1p(-x);
		const double d81 = 1.0/81.0;

		const double nf = static_cast<double>(_nf);

		double res = 1641.1 - 3135.0*x + 243.6*std::pow(x, 2.0) - 522.1*std::pow(x, 3.0)
			+ 128.*d81*std::pow(dl, 4.0) + 2400.*d81*std::pow(dl, 3.0)
			+ 294.9*std::pow(dl, 2.0) + 1258.0*dl
			+ 714.1*dlm + dl*dlm*(563.9 + 256.8*dl)
			+ nf * ( -197.0 + 381.1*x + 72.94*std::pow(x, 2.0) + 44.79*std::pow(x, 3.0)
					 - 192.0*d81*std::pow(dl, 3.0) - 2608.0*d81*std::pow(dl, 2.0) - 152.6*dl
					 - 5120.0*d81*dlm - 56.66*dl*dlm - 1.497*x*std::pow(dl, 3.0) )
			+ std::pow(nf, 2.0)*( 32.0*x*dl/(1.0-x) * (3.0*dl + 10.0) + 64.0
								  + (48.0*std::pow(dl, 2.0) + 352.0*dl + 384.0)*(1.0-x) ) * d81;

		return res/8.0;
	}
	double P2nsp::calcPlus(double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		
		double res = (1174.898 - nf*183.187 - pow(nf, 2)*(64.0/81.0));
		
		return res/8.0;
	}
	double P2nsp::calcDelta(double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		// const double dl1 = std::log1p(-x);

		double res = 1295.624 - 0.24 - nf*(173.938 - 0.011) + std::pow(nf, 2)*1.13067;
		
		return res/8.0;
	}


	double P2nsm::calcRegular(double x) const
	{
		const double dl = std::log(x);
		const double dl2 = dl*dl;
		const double dl3 = dl2*dl;
		const double dl4 = dl3*dl;
		const double dlm = std::log1p(-x);
		const double d81 = 1.0/81.0;
		const double x2 = x*x;
		const double x3 = x2*x;

		const double nf = static_cast<double>(_nf);

		double res = 1860.2 - 3505.0*x + 297.0*x2 - 433.2*x3
			+ 116.0*d81*dl4 + 2880.0*d81*dl3
			+ 399.2*dl2 + 1465.2*dl
			+ 714.1*dlm + dl*dlm*(684.0 + 251.2*dl)
			+ nf*( -216.62 + 406.5*x + 77.89*x2 + 34.76*x3
				   - 256.0*d81*dl3 - 3216.0*d81*dl2 - 172.69*dl 
				   - 5120.*d81*dlm - 65.43*dl*dlm - 1.136*x*dl3 )
			+ std::pow(nf, 2.0)*( 32.0*x*dl/(1.0-x) * (3.0*dl + 10.0) + 64.0
								  + (48.0*dl2 + 352.0*dl + 384.0)*(1.0-x) ) * d81;

		return res/8.0;
	}
	double P2nsm::calcPlus(double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		
		double res = (1174.898 - nf*183.187 - pow(nf, 2)*(64.0/81.0));
		
		return res/8.0;
	}
	double P2nsm::calcDelta(double x) const
	{
		UNUSED(x);

		const double nf = static_cast<double>(_nf);

		double res = (1295.624 - 0.154) - nf*(173.938 - 0.005) + std::pow(nf, 2.0)*(1.13067);

		return res/8.0;
	}


	double P2nsv::calcRegular(double x) const
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
	double P2nsv::calcPlus(double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		
		double res = (1174.898 - nf*183.187 - pow(nf, 2)*(64.0/81.0));
		
		return res/8.0;
	}
	double P2nsv::calcDelta(double x) const
	{
		UNUSED(x);

		const double nf = static_cast<double>(_nf);

		double res = (1295.624 - 0.154) - nf*(173.938 - 0.005) + std::pow(nf, 2.0)*(1.13067);

		return res/8.0;
	}

	double P2ps::calcRegular(double x) const
	{
		double dl  = std::log(x);
        double dl1 = std::log1p(-x);

		double NF = static_cast<double>(_nf);

		
		double resa = - 3584./(27.*x) * dl - 506.0/ x + 160./27. * std::pow(dl, 4)
			- 400./9. * std::pow(dl, 3) + 131.4 * std::pow(dl, 2) - 661.6 * dl
			- 5.926  * std::pow(dl1, 3) - 9.751 * std::pow(dl1, 2) - 72.11 * dl1
			+ 177.4 + 392.9 * x - 101.4 * std::pow(x, 2) - 57.04 * dl*dl1;
		double resb =  256./(81.*x) + 32./27. * std::pow(dl, 3) + 17.89 * std::pow(dl, 2)
			+ 61.75 * dl + 1.778 * std::pow(dl1, 2) + 5.944 * dl1 + 100.1
			- 125.2 * x + 49.26 * std::pow(x, 2) - 12.59 * std::pow(x, 3) 
			- 1.889 * dl*dl1;
			
		double res = (1.0-x)*NF*(resa + NF*resb);
		return res/8.0;
	}

	double P2qq::calcRegular(double x) const
	{
		double dl  = std::log(x);
        double dl1 = std::log1p(-x);
		double d81 = 1.0/81.0;

		double NF = static_cast<double>(_nf);


		double res1 =   1641.1 - 3135.* x + 243.6 * std::pow(x, 2) - 522.1 * std::pow(x, 3)
                 + 128.*d81 * std::pow(dl, 4) + 2400.*d81 * std::pow(dl, 3)
                 + 294.9 * std::pow(dl, 2) + 1258.* dl
                 + 714.1 * dl1 + dl*dl1 * (563.9 + 256.8 * dl)
             + NF * ( -197.0 + 381.1 * x + 72.94 * std::pow(x, 2) + 44.79 * std::pow(x, 3)
                 - 192.*d81 * std::pow(dl, 3)  - 2608.*d81 * std::pow(dl, 2) - 152.6 * dl
                 - 5120.*d81 * dl1 - 56.66 * dl*dl1 - 1.497 * x*std::pow(dl, 3) )
			+ std::pow(NF, 2) * ( 32.* x*dl/(1.-x) * (3.* dl + 10.) + 64.
						 + (48.* std::pow(dl, 2) + 352.* dl + 384.) * (1.-x) ) * d81;

		
		double res2a = - 3584./(27.*x) * dl - 506.0/ x + 160./27. * std::pow(dl, 4)
			- 400./9. * std::pow(dl, 3) + 131.4 * std::pow(dl, 2) - 661.6 * dl
			- 5.926  * std::pow(dl1, 3) - 9.751 * std::pow(dl1, 2) - 72.11 * dl1
			+ 177.4 + 392.9 * x - 101.4 * std::pow(x, 2) - 57.04 * dl*dl1;
		double res2b =  256./(81.*x) + 32./27. * std::pow(dl, 3) + 17.89 * std::pow(dl, 2)
			+ 61.75 * dl + 1.778 * std::pow(dl1, 2) + 5.944 * dl1 + 100.1
			- 125.2 * x + 49.26 * std::pow(x, 2) - 12.59 * std::pow(x, 3) 
			- 1.889 * dl*dl1;
			
		double res2 = (1.0-x)*NF*(res2a + NF*res2b);

		return (res1 + res2)/8.0;
		
	}
	double P2qq::calcPlus(double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		
		double res = (1174.898 - nf*183.187 - std::pow(nf, 2)*(64.0/81.0));
		
		return res/8.0;
	}
	double P2qq::calcDelta(double x) const
	{
		UNUSED(x);
		
		const double nf = static_cast<double>(_nf);
		// const double dl1 = std::log1p(-x);

		double res = 1295.624 - 0.24 - nf*(173.938 - 0.011) + std::pow(nf, 2)*1.13067;
		
		return res/8.0;
	}


	double P2qg::calcRegular(double x) const
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


	double P2gq::calcRegular(double x) const
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


	double P2gg::calcRegular(double x) const
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
	double P2gg::calcPlus(double x) const
	{
		UNUSED(x);
		
		const double Nf = static_cast<double>(_nf);
		
		double res = 2643.521 - Nf*412.172 - Nf*Nf*(16.0/9.0);

		return res/8.0;
	}
	double P2gg::calcDelta(double x) const
	{
		UNUSED(x);
		
		const double Nf = static_cast<double>(_nf);

		//const double dl1 = std::log1p(-x);

		double res = (4425.448 + 0.446) - Nf*(528.720 + 0.003) + Nf*Nf*(6.4630);
		
		return res/8.0;
	}




	
	double P3nsp::calcRegular(double x) const
	{	
		const double x2   = x*x;
		const double x3   = x2*x;
		const double omx  = 1.0-x;
		const double dm   = 1.0/omx;
		const double dl   = std::log(x);
		const double dl2  = dl*dl;
		const double dl3  = dl2*dl;
		const double dl4  = dl3*dl;
		const double dl5  = dl4*dl;
		const double dl6  = dl5*dl;
		const double dlm  = std::log1p(-x);
		const double dlm2 = dlm*dlm;
		const double dlm3 = dlm2*dlm;

		// Leading large-n_c, nf^0 and nf^1, parametrized
		const double p3nsa0  =
			2.5e4*(omx*(3.5254 + 8.6935*x - 1.5051*x2 + 1.8300*x3 )
					   + 11.883*x*dl - 0.09066*x*dl2 + 11.410*omx*dlm + 13.376*dl*dlm)
			+ 5.167133e+4*dl + 1.712095e+4*dl2 + 2.863226e+3*dl3 + 2.978255e+2*dl4
			+ 1.6e+1 * dl5 + 5.e-1 * dl6 - 2.973385e+4 + 1.906980e+4 * dlm;
		const double p3nsa1  =
			2.5e+4 * ( omx * ( - 0.74077 + 1.4860 * x - 0.23631 * x2 + 0.31584 * x3 )
					   + 2.5251 * omx * dlm + 2.5203 * dl * dlm + 2.2242 * x * dl
					   - 0.02460 * x * dl2 + 0.00310 * x * dl3 )
			- 9.239374e+3 * dl - 2.917312e+3 * dl2 - 4.305308e+2 *dl3 - 3.6e+1 * dl4
			- 4.0/3.0 * dl5 + 8.115605e+3 - 3.079761e+3 * dlm;

		// Nonleading large-n_c, nf^0 and nf^1: two approximations
		const double p3npa01 =
			3948.16 * omx - 2464.61 * ( 2.0*x - x2 ) * omx - 1839.44 * dl2 - 402.156 * dl3
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
			- ( 2.633745e-1 * dm - 1.31687e-1 * ( 1.0 + x ) ) * dl3;

		// Assembly
		double nf = static_cast<double>(_nf);
		const double p3nspai = p3nsa0 + nf * p3nsa1 + nf * nf * p3nspa2 + nf * nf * nf * p3nsa3;
	    double res = std::numeric_limits<double>::max();
		if (_imod == 1)
			res = p3nspai + p3npa01 + nf * p3npa11;
		else if (_imod == 2)
			res = p3nspai + p3npa02 + nf * p3npa12;
		else
			res = p3nspai + 0.5 * ( ( p3npa01 + p3npa02 ) + nf * ( p3npa11 + p3npa12 ) );

		return res/16.0;
	}
	double P3nsp::calcPlus(double x) const
	{
		UNUSED(x);
		
		double Nf = static_cast<double>(_nf);

		const double a4qi =
			2.120902e+4
		  - 5.179372e+3*Nf
		  + 1.955772e+2*Nf*Nf
		  + 3.272344e+0*Nf*Nf*Nf;
		const double a4ap1 = - 507.152 + 7.33927*Nf;
		const double a4ap2 = - 505.209 + 7.53662*Nf;

	    double res = std::numeric_limits<double>::max();
		if (_imod == 1)
			res = a4qi + a4ap1;
		else if (_imod == 2)
			res = a4qi + a4ap2;
		else
			res = a4qi + 0.5*(a4ap1+a4ap2);

		return res/16.0;
	}
	double P3nsp::calcDelta(double x) const
	{
		UNUSED(x);
		
		double Nf = static_cast<double>(_nf);

		// this is for the coefficients of a log(1-x) component,
		// which is present in the event one does some antiderivatives first
		// related to the plus distribution
		// --> anti-derivative of 1/(1-x) is log(1-x)
		// those pieces end up here somehow, but that is now how we do it in candia
		// the pieces are kept here anyways just in case
		/*
		const double a4qi  =
			2.120902e+4
			- 5.179372e+3 * _nf
			+ 1.955772e+2 * _nf * _nf
			+ 3.272344e+0 * _nf * _nf * _nf;
		const double a4ap1 = - 507.152 + 7.33927 * _nf;
		const double a4ap2 = - 505.209 + 7.53662 * _nf;
		*/

		const double b4qi =
			2.579609e+4 + 0.08
		  - (5.818637e+3+0.97)   *Nf
		  + (1.938554e+2+0.0037) *Nf*Nf
		  +  3.014982e+0         *Nf*Nf*Nf;
		const double b4ap1 = - 2405.03 + 267.965 * Nf;
		const double b4ap2 = - 2394.47 + 269.028 * Nf;

	    double res = std::numeric_limits<double>::max();
		if (_imod == 1)
			res = b4qi + b4ap1;
		else if (_imod == 2)
			res = b4qi + b4ap2;
		else
			res = b4qi + 0.5 * ( b4ap1 + b4ap2 );

		return res/16.0;
	}


	double P3nsm::calcRegular(double x) const
	{
		const double x2   = x*x;
		const double x3   = x2*x;
		const double omx  = 1.0-x;
		const double dm   = 1.0/omx;
		const double dl   = std::log(x);
		const double dl2  = dl*dl;
		const double dl3  = dl2*dl;
		const double dl4  = dl3*dl;
		const double dl5  = dl4*dl;
		const double dl6  = dl5*dl;
		const double dlm  = std::log1p(-x);
		const double dlm2 = dlm*dlm;
		const double dlm3 = dlm2*dlm;

		// Leading large-n_c, nf^0 and nf^1, parametrized
		const double p3nsa0  =
			2.5e+4 * ( omx * ( 3.5254 + 8.6935 * x - 1.5051 * x2 + 1.8300 * x3 )
					   + 11.883 * x * dl - 0.09066 * x * dl2 + 11.410 * omx * dlm + 13.376  * dl * dlm )
			+ 5.167133e+4 * dl + 1.712095e+4 * dl2 + 2.863226e+3 * dl3 + 2.978255e+2 * dl4
			+ 1.6e+1 * dl5 + 5.0e-1 * dl6 - 2.973385e+4 + 1.906980e+4 * dlm;

		const double p3nsa1  =
			2.5e+4 * ( omx * ( - 0.74077 + 1.4860 * x - 0.23631 * x2 + 0.31584 * x3 )
					   + 2.5251 * omx * dlm + 2.5203 * dl * dlm + 2.2242 * x * dl
					   - 0.02460 * x * dl2 + 0.00310 * x * dl3 )
			- 9.239374e+3 * dl - 2.917312e+3 * dl2 - 4.305308e+2 *dl3 - 3.6e+1 * dl4
			- 4. / 3. * dl5 + 8.115605e+3 - 3.079761e+3 * dlm;

		// Nonleading large-n_c, nf^0 and nf^1: two approximations
		const double p3nma01 =
			( 5992.88 * ( 1.0 + 2.0 * x ) + 31321.44 * x2 ) * omx + 511.228 - 1618.07 * dl + 2.25480 * dl3
			+ 31897.82 * dlm * omx + 4653.76 * dlm2 * omx + 4.964335e-1 * ( dl6 + 6.0 * dl5 )
			- 2.601749e+3 - 2.118867e+3 * dlm;
		const double p3nma02 =
			( 4043.59 - 15386.6 * x ) * x * omx + 502.481 + 1532.96  * dl2 + 31.6023 * dl3
			- 3997.39  * dlm * omx + 511.567 * dlm3 * omx + 4.964335e-1 * ( dl6 + 18.0 * dl5 )
			- 2.601749e+3 - 2.118867e+3 * dlm;

		const double p3nma11 =
			( 114.457 * ( 1.0 + 2.0 * x ) + 2570.73 * x2 ) * omx - 7.08645 - 127.012 * dl2 + 2.69618 * dl4
			+ 1856.63 * dlm * omx + 440.17 * dlm2 * omx + 3.121643e+2 + 3.379310e+2 * dlm;
		const double p3nma12 =
			( - 335.995 * ( 2.0 + x ) - 1605.91 * x2 ) * omx - 7.82077 - 9.76627 * dl2 + 0.14218 * dl5
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
			- ( 1.316872 * ( dm + 1.0e-1) - 1.448560 * x ) * dl2
			- ( 2.633744e-1 * dm - 1.31687e-1 * ( 1.0 + x ) ) * dl3;

		// Assembly
		const double Nf = static_cast<double>(_nf);
		const double p3nsmai = p3nsa0 + Nf * p3nsa1 + Nf * Nf * p3nsma2 + Nf * Nf * Nf * p3nsa3;
	    double res = std::numeric_limits<double>::max();
		if (_imod == 1)
			res = p3nsmai + p3nma01 + Nf * p3nma11;
		else if (_imod == 2)
			res = p3nsmai + p3nma02 + Nf * p3nma12;
		else
			res = p3nsmai + 0.5*((p3nma01 + p3nma02) + Nf*(p3nma11 + p3nma12));

		return res/16.0;
	}
	double P3nsm::calcPlus(double x) const
	{
		UNUSED(x);
		double Nf = static_cast<double>(_nf);

		const double a4qi  =
			2.120902e+4
			- 5.179372e+3 * Nf
			+ 1.955772e+2 * Nf * Nf
			+ 3.272344e+0 * Nf * Nf * Nf;
		const double a4ap1 = - 511.228 + 7.08645 * Nf;
		const double a4ap2 = - 502.481 + 7.82077 * Nf;

		double res = std::numeric_limits<double>::max();
		if (_imod == 1)
			res = a4qi + a4ap1;
		else if (_imod == 2)
			res = a4qi + a4ap2;
		else
			res = a4qi + 0.5*(a4ap1+a4ap2);

		return res/16.0;
	}
	double P3nsm::calcDelta(double x) const
	{
		UNUSED(x);
		double Nf = static_cast<double>(_nf);

		/*
		const double a4qi  =
			2.120902e+4
			- 5.179372e+3 * _nf
			+ 1.955772e+2 * _nf * _nf
			+ 3.272344e+0 * _nf * _nf * _nf;
		const double a4ap1 = - 511.228 + 7.08645 * _nf;
		const double a4ap2 = - 502.481 + 7.82077 * _nf;
		*/

		const double b4qi =
			2.579609e+4 + 0.08
			- ( 5.818637e+3 + 0.97 )   * Nf
			+ ( 1.938554e+2 + 0.0037 ) * Nf * Nf
			+   3.014982e+0 * Nf * Nf * Nf;
		const double b4ap1 = - 2426.05  + 266.674 * Nf - 0.05 * Nf;
		const double b4ap2 = - 2380.255 + 270.518 * Nf - 0.05 * Nf;

		double res = std::numeric_limits<double>::max();
		if (_imod == 1)
			res = b4qi + b4ap1;
		else if (_imod == 2)
			res = b4qi + b4ap2;
		else
			res = b4qi + 0.5*(b4ap1+b4ap2);

		return res/16.0;
	}

	double P3nsv::calcRegular(double x) const
	{
		const double Nf = static_cast<double>(_nf);
		double res1 = std::numeric_limits<double>::max();
		double res2 = std::numeric_limits<double>::max();
		{
			const double x2   = x * x;
			const double omx  = 1.0 - x;
			const double dl   = std::log(x);
			const double dl2  = dl * dl;
			const double dl3  = dl2 * dl;
			const double dl4  = dl3 * dl;
			const double dl5  = dl4 * dl;
			const double dl6  = dl5 * dl;
			const double dlm  = std::log(omx);
			const double dlm2 = dlm * dlm;
			const double dlm3 = dlm2 * dlm;

			// nf^1: two approximations
			const double p3nsa11 =
				omx * x * ( 4989.2 - 1607.73 * x ) + 3687.6 * dl + 3296.6 * dl2 + 1271.11* dl3
				+ 533.44 * dl4 + 97.27 *  dl5 + 4 * dl6 + 60.40 * omx * dlm2 + 4.685 * omx * dlm3;
			const double p3nsa12 =
				1030.79 * omx * x + 1266.77 * omx * ( 2.0 - x2 ) + 2987.83 * dl + 273.05 * dl2 - 923.48 * dl3
				- 236.76 * dl4 - 33.886 * dl5 - 4.0 * dl6 - 254.63 * omx * dlm - 0.28953 * omx * dlm3;

			// nf^2 (parametrized)
			const double p3nssa2 =
				2.5e+2 * ( omx * ( - 4.7656 + 1.6908 * x + 0.1703 * x2 )
						- 0.41652 * x *dl + 0.90777 * x * dl2 + 0.12478 * x * dl3
						+ 0.17155 * omx * dlm + 0.17191  * dl * dlm )
				- 6.473971e+2 * dl - 6.641219e+1 * dl2 - 5.353347 * dl3 - 5.925926 * dl4
				- 3.950617e-1 * dl5 + 1.970002e+1 * omx * dlm - 3.435474 * omx * dlm2;

			if (_imod == 1)
				res1 = Nf * p3nsa11 + Nf * Nf * p3nssa2;
			else if (_imod == 2)
				res1 = Nf * p3nsa12 + Nf * Nf * p3nssa2;
			else
				res1 = 0.5 *Nf * ( p3nsa11 + p3nsa12 ) + Nf * Nf * p3nssa2;
		}

		{
			const double x2   = x * x;
			const double x3   = x2 * x;
			const double omx  = 1.0 - x;
			const double dm   = 1.0 / omx;
			const double dl   = std::log(x);
			const double dl2  = dl * dl;
			const double dl3  = dl2 * dl;
			const double dl4  = dl3 * dl;
			const double dl5  = dl4 * dl;
			const double dl6  = dl5 * dl;
			const double dlm  = std::log(omx);
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
				( 5992.88 * ( 1. + 2. * x ) + 31321.44 * x2 ) * omx + 511.228 - 1618.07 * dl + 2.25480 * dl3
				+ 31897.82 * dlm * omx + 4653.76 * dlm2 * omx + 4.964335e-1 * ( dl6 + 6. * dl5 )
				- 2.601749e+3 - 2.118867e+3 * dlm;
			const double p3nma02 =
				( 4043.59 - 15386.6 * x ) * x * omx + 502.481 + 1532.96  * dl2 + 31.6023 * dl3
				- 3997.39  * dlm * omx + 511.567 * dlm3 * omx + 4.964335e-1 * ( dl6 + 18. * dl5 )
				- 2.601749e+3 - 2.118867e+3 * dlm;

			const double p3nma11 =
				( 114.457 * ( 1. + 2. * x ) + 2570.73 * x2 ) * omx - 7.08645 - 127.012 * dl2 + 2.69618 * dl4
				+ 1856.63 * dlm * omx + 440.17 * dlm2 * omx + 3.121643e+2 + 3.379310e+2 * dlm;
			const double p3nma12 =
				( - 335.995 * ( 2. + x ) - 1605.91 * x2 ) * omx - 7.82077 - 9.76627 * dl2 + 0.14218 * dl5
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
				- ( 1.316872 * ( dm + 1.e-1) - 1.448560 * x ) * dl2
				- ( 2.633744e-1 * dm - 1.31687e-1 * ( 1. + x ) ) * dl3;

			// Assembly
			const double p3nsmai = p3nsa0 + Nf * p3nsa1 + Nf * Nf * p3nsma2 + Nf * Nf * Nf * p3nsa3;			
			if (_imod == 1)
				res2 = p3nsmai + p3nma01 + Nf * p3nma11;
			else if (_imod == 2)
				res2 = p3nsmai + p3nma02 + Nf * p3nma12;
			else
				res2 = p3nsmai + 0.5 * ( ( p3nma01 + p3nma02 ) + Nf * ( p3nma11 + p3nma12 ) );
		}

		return (res1+res2)/16.0;
	}
	double P3nsv::calcPlus(double x) const
	{
		UNUSED(x);
		double Nf = static_cast<double>(_nf);

		const double a4qi  = 2.120902e+4
			- 5.179372e+3*Nf
			+ 1.955772e+2*Nf*Nf
			+ 3.272344e+0*Nf*Nf*Nf;
		const double a4ap1 = -511.228 + 7.08645*Nf;
		const double a4ap2 = -502.481 + 7.82077*Nf;

	    double res = std::numeric_limits<double>::max();
		if (_imod == 1)
			res = a4qi + a4ap1;
		else if (_imod == 2)
			res = a4qi + a4ap2;
		else
			res = a4qi + 0.5*(a4ap1+a4ap2);

		return res/16.0;
	}
	double P3nsv::calcDelta(double x) const
	{
		UNUSED(x);
		double Nf = static_cast<double>(_nf);

		/*
		const double a4qi  =
			2.120902e+4
			- 5.179372e+3 * _nf
			+ 1.955772e+2 * _nf * _nf
			+ 3.272344e+0 * _nf * _nf * _nf;
		const double a4ap1 = - 511.228 + 7.08645 * _nf;
		const double a4ap2 = - 502.481 + 7.82077 * _nf;
		*/

		const double b4qi =
			2.579609e+4 + 0.08
			- ( 5.818637e+3 + 0.97 )   * Nf
			+ ( 1.938554e+2 + 0.0037 ) * Nf * Nf
			+   3.014982e+0 * Nf * Nf * Nf;
		const double b4ap1 = - 2426.05  + 266.674 * Nf - 0.05 * Nf;
		const double b4ap2 = - 2380.255 + 270.518 * Nf - 0.05 * Nf;

	    double res = std::numeric_limits<double>::max();
		if (_imod == 1)
			res = b4qi + b4ap1;
		else if (_imod == 2)
			res = b4qi + b4ap2;
		else
			res = b4qi + 0.5*(b4ap1+b4ap2);

		return res/16.0;
	}




	double P3ps::calcRegular(double x) const
	{
		const double Nf = static_cast<double>(_nf);
		const double Nf2     = Nf*Nf;
		const double Nf3     = Nf*Nf2;
		const double xm   = 1.0 / x;
		const double x1   = 1.0 - x;
		const double dl   = std::log(x);
		const double dl2  = dl * dl;
		const double dl3  = dl * dl2;
		const double dl4  = dl * dl3;
		const double dl5  = dl * dl4;
		const double dl6  = dl * dl5;
		const double dlm  = std::log1p(-x);
		const double dlm2 = dlm * dlm;
		const double dlm3 = dlm * dlm2;
		const double dlm4 = dlm * dlm3;

		// Known large-x coefficients
		const double x1L4cff = - 5.6460905e1 * Nf + 3.6213992   * Nf2;
		const double x1L3cff = - 2.4755054e2 * Nf + 4.0559671e1 * Nf2 - 1.5802469 * Nf3;
		const double y1L4cff = - 1.3168724e1 * Nf;
		const double y1L3cff = - 1.9911111e2 * Nf + 1.3695473e1 * Nf2;

		// Known small-x coefficients
		const double bfkl1   =   1.7492273e3 * Nf;
		const double x0L6cff = - 7.5061728   * Nf + 7.9012346e-1 * Nf2;
		const double x0L5cff =   2.8549794e1 * Nf + 3.7925926    * Nf2;
		const double x0L4cff = - 8.5480010e2 * Nf + 7.7366255e1  * Nf2 - 1.9753086e-1 * Nf3;

		// The resulting part of the function
		const double P3ps01 =
			+ bfkl1 * dl2 * xm
			+ x0L6cff * dl6
			+ x0L5cff * dl5
			+ x0L4cff * dl4
			+ x1L3cff * x1 * dlm3
			+ x1L4cff * x1 * dlm4
			+ y1L3cff * x1 * x1 * dlm3
			+ y1L4cff * x1 * x1 * dlm4;

		// The selected approximations for nf = 3, 4, 5
		double P3psApp1 = P3ps01;
		double P3psApp2 = P3ps01;
		if (_nf <= 3) {
			P3psApp1 +=
				+ 67731.  * x1 * dl * xm
				+ 274100. * x1 * xm
				- 104493. * x1 * ( 1. + 2. * x )
				+ 34403.  * x1 * x * x
				+ 353656. * x1 * dl
				+ 10620.  * dl2
				+ 40006.  * dl3
				- 7412.1  * x1 * dlm
				- 2365.1  * x1 * dlm2
				+ 1533.0  * x1 * x1 * dlm2;
			P3psApp2 +=
				+ 54593.  * x1 * dl * xm
				+ 179748. * x1 * xm
				- 195263. * x1
				+ 12789.  * x1 * x * ( 1. + x )
				+ 4700.0  * x1 * dl
				- 103604. * dl2
				- 2758.3  * dl3
				- 2801.2  * x1 * dlm
				- 1986.9  * x1 * dlm2
				- 6005.9  * x1 * x1 * dlm2;
		} else if (_nf == 4) {
			P3psApp1 +=
				+ 90154.  * x1 * dl *xm
				+ 359084. * x1 * xm
				- 136319. * x1 * ( 1. + 2. * x )
				+ 45379.  * x1 * x * x
				+ 461167. * x1 * dl
				+ 13869.  * dl2
				+ 52525.  * dl3
				- 7498.2  * x1 * dlm
				- 2491.5  * x1 * dlm2
				+ 1727.2  * x1 * x1 * dlm2;
			P3psApp2 +=
				+ 72987.  * x1 * dl * xm
				+ 235802. * x1 * xm
				- 254921. * x1
				+ 17138.  * x1 * x * ( 1. + x )
				+ 5212.9  * x1 * dl
				- 135378. * dl2
				- 3350.9  * dl3
				- 1472.7  * x1 * dlm
				- 1997.2  * x1 * dlm2
				- 8123.3  * x1 * x1 * dlm2;
		} else if (_nf >= 5) {
			P3psApp1 +=
				+ 112481. * x1 * dl * xm
				+ 440555. * x1 * xm
				- 166581. * x1 * ( 1. + 2. * x )
				+ 56087.  * x1 * x * x
				+ 562992. * x1 * dl
				+ 16882.  * dl2
				+ 64577.  * dl3
				- 6570.1  * x1 * dlm
				- 2365.7  * x1 * dlm2
				+ 1761.7  * x1 * x1 * dlm2;
			P3psApp2 +=
				+ 91468.  * x1 * dl * xm
				+ 289658. * x1 * xm
				- 311749. * x1
				+ 21521.  * x1 * x * ( 1. + x )
				+ 4908.9 * x1 * dl
				- 165795. * dl2
				- 3814.9 * dl3
				+ 804.5 * x1 * dlm
				- 1760.8 * x1 * dlm2
				- 10295.  * x1 * x1 * dlm2;
		}

		// We return (for now) one of the two error-band boundaries or the
		// present best estimate, their average
		double res = std::numeric_limits<double>::max()/16.0;
		if (_imod == 1)
			res = P3psApp1;
		else if (_imod == 2)
			res = P3psApp2;
		else
			res = 0.5 * ( P3psApp1 + P3psApp2 );

		return res/16.0;
	}

	double P3qq::calcRegular(double x) const
	{
		const double Nf = static_cast<double>(_nf);
		double res1 = std::numeric_limits<double>::max(),
			   res2 = std::numeric_limits<double>::max();

		// P3nsp
		{
			const double x2   = x * x;
			const double x3   = x2 * x;
			const double omx  = 1.0 - x;
			const double dm   = 1.0 / omx;
			const double dl   = std::log(x);
			const double dl2  = dl * dl;
			const double dl3  = dl2 * dl;
			const double dl4  = dl3 * dl;
			const double dl5  = dl4 * dl;
			const double dl6  = dl5 * dl;
			const double dlm  = std::log1p(-x);
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
				- 4.0/3.0 * dl5 + 8.115605e+3 - 3.079761e+3 * dlm;

			// Nonleading large-n_c, nf^0 and nf^1: two approximations
			const double p3npa01 =
				3948.16 * omx - 2464.61 * ( 2.0*x - x2 ) * omx - 1839.44 * dl2 - 402.156 * dl3
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
				- ( 1.316872 * ( dm + 1.e-1) - 1.448560 * x ) * dl2
				- ( 2.633745e-1 * dm - 1.31687e-1 * ( 1. + x ) ) * dl3;

			// Assembly
			const double p3nspai = p3nsa0 + Nf * p3nsa1 + Nf * Nf * p3nspa2 + Nf * Nf * Nf * p3nsa3;
			if (_imod == 1)
				res1 = p3nspai + p3npa01 + Nf * p3npa11;
			else if (_imod == 2)
				res1 = p3nspai + p3npa02 + Nf * p3npa12;
			else
				res1 = p3nspai + 0.5 * ( ( p3npa01 + p3npa02 ) + Nf * ( p3npa11 + p3npa12 ) );
		}

		// P3ps
		{
			const double Nf2     = Nf*Nf;
			const double Nf3     = Nf*Nf2;
			const double xm   = 1.0 / x;
			const double x1   = 1.0 - x;
			const double dl   = std::log(x);
			const double dl2  = dl * dl;
			const double dl3  = dl * dl2;
			const double dl4  = dl * dl3;
			const double dl5  = dl * dl4;
			const double dl6  = dl * dl5;
			const double dlm  = std::log1p(-x);
			const double dlm2 = dlm * dlm;
			const double dlm3 = dlm * dlm2;
			const double dlm4 = dlm * dlm3;

			// Known large-x coefficients
			const double x1L4cff = - 5.6460905e1 * Nf + 3.6213992   * Nf2;
			const double x1L3cff = - 2.4755054e2 * Nf + 4.0559671e1 * Nf2 - 1.5802469 * Nf3;
			const double y1L4cff = - 1.3168724e1 * Nf;
			const double y1L3cff = - 1.9911111e2 * Nf + 1.3695473e1 * Nf2;

			// Known small-x coefficients
			const double bfkl1   =   1.7492273e3 * Nf;
			const double x0L6cff = - 7.5061728   * Nf + 7.9012346e-1 * Nf2;
			const double x0L5cff =   2.8549794e1 * Nf + 3.7925926    * Nf2;
			const double x0L4cff = - 8.5480010e2 * Nf + 7.7366255e1  * Nf2 - 1.9753086e-1 * Nf3;

			// The resulting part of the function
			const double P3ps01 =
				+ bfkl1 * dl2 * xm
				+ x0L6cff * dl6
				+ x0L5cff * dl5
				+ x0L4cff * dl4
				+ x1L3cff * x1 * dlm3
				+ x1L4cff * x1 * dlm4
				+ y1L3cff * x1 * x1 * dlm3
				+ y1L4cff * x1 * x1 * dlm4;

			// The selected approximations for nf = 3, 4, 5
			double P3psApp1 = P3ps01;
			double P3psApp2 = P3ps01;
			if (_nf <= 3) {
				P3psApp1 +=
					+ 67731.  * x1 * dl * xm
					+ 274100. * x1 * xm
					- 104493. * x1 * ( 1. + 2. * x )
					+ 34403.  * x1 * x * x
					+ 353656. * x1 * dl
					+ 10620.  * dl2
					+ 40006.  * dl3
					- 7412.1  * x1 * dlm
					- 2365.1  * x1 * dlm2
					+ 1533.0  * x1 * x1 * dlm2;
				P3psApp2 +=
					+ 54593.  * x1 * dl * xm
					+ 179748. * x1 * xm
					- 195263. * x1
					+ 12789.  * x1 * x * ( 1. + x )
					+ 4700.0  * x1 * dl
					- 103604. * dl2
					- 2758.3  * dl3
					- 2801.2  * x1 * dlm
					- 1986.9  * x1 * dlm2
					- 6005.9  * x1 * x1 * dlm2;
			} else if (_nf == 4) {
				P3psApp1 +=
					+ 90154.  * x1 * dl *xm
					+ 359084. * x1 * xm
					- 136319. * x1 * ( 1. + 2. * x )
					+ 45379.  * x1 * x * x
					+ 461167. * x1 * dl
					+ 13869.  * dl2
					+ 52525.  * dl3
					- 7498.2  * x1 * dlm
					- 2491.5  * x1 * dlm2
					+ 1727.2  * x1 * x1 * dlm2;
				P3psApp2 +=
					+ 72987.  * x1 * dl * xm
					+ 235802. * x1 * xm
					- 254921. * x1
					+ 17138.  * x1 * x * ( 1. + x )
					+ 5212.9  * x1 * dl
					- 135378. * dl2
					- 3350.9  * dl3
					- 1472.7  * x1 * dlm
					- 1997.2  * x1 * dlm2
					- 8123.3  * x1 * x1 * dlm2;
			} else if (_nf >= 5) {
				P3psApp1 +=
					+ 112481. * x1 * dl * xm
					+ 440555. * x1 * xm
					- 166581. * x1 * ( 1. + 2. * x )
					+ 56087.  * x1 * x * x
					+ 562992. * x1 * dl
					+ 16882.  * dl2
					+ 64577.  * dl3
					- 6570.1  * x1 * dlm
					- 2365.7  * x1 * dlm2
					+ 1761.7  * x1 * x1 * dlm2;
				P3psApp2 +=
					+ 91468.  * x1 * dl * xm
					+ 289658. * x1 * xm
					- 311749. * x1
					+ 21521.  * x1 * x * ( 1. + x )
					+ 4908.9 * x1 * dl
					- 165795. * dl2
					- 3814.9 * dl3
					+ 804.5 * x1 * dlm
					- 1760.8 * x1 * dlm2
					- 10295.  * x1 * x1 * dlm2;
			}

			// We return (for now) one of the two error-band boundaries or the
			// present best estimate, their average
			if (_imod == 1)
				res2 = P3psApp1/16.0;
			else if (_imod == 2)
				res2 = P3psApp2;
			else
				res2 = 0.5 * ( P3psApp1 + P3psApp2 );
		}

		return (res1+res2)/16.0;
	}
	double P3qq::calcPlus(double x) const
	{
		UNUSED(x);
		
		double Nf = static_cast<double>(_nf);

		const double a4qi =
			2.120902e+4
		  - 5.179372e+3*Nf
		  + 1.955772e+2*Nf*Nf
		  + 3.272344e+0*Nf*Nf*Nf;
		const double a4ap1 = - 507.152 + 7.33927*Nf;
		const double a4ap2 = - 505.209 + 7.53662*Nf;

	    double res = std::numeric_limits<double>::max();
		if (_imod == 1)
			res = a4qi + a4ap1;
		else if (_imod == 2)
			res = a4qi + a4ap2;
		else
			res = a4qi + 0.5*(a4ap1+a4ap2);

		return res/16.0;
	}
	double P3qq::calcDelta(double x) const
	{
		UNUSED(x);
		
		double Nf = static_cast<double>(_nf);

		// this is for the coefficients of a log(1-x) component,
		// which is present in the event one does some antiderivatives first
		// related to the plus distribution
		// --> anti-derivative of 1/(1-x) is log(1-x)
		// those pieces end up here somehow, but that is now how we do it in candia
		// the pieces are kept here anyways just in case
		/*
		const double a4qi  =
			2.120902e+4
			- 5.179372e+3 * _nf
			+ 1.955772e+2 * _nf * _nf
			+ 3.272344e+0 * _nf * _nf * _nf;
		const double a4ap1 = - 507.152 + 7.33927 * _nf;
		const double a4ap2 = - 505.209 + 7.53662 * _nf;
		*/

		const double b4qi =
			2.579609e+4 + 0.08
		  - (5.818637e+3+0.97)   *Nf
		  + (1.938554e+2+0.0037) *Nf*Nf
		  +  3.014982e+0         *Nf*Nf*Nf;
		const double b4ap1 = - 2405.03 + 267.965 * Nf;
		const double b4ap2 = - 2394.47 + 269.028 * Nf;

	    double res = std::numeric_limits<double>::max();
		if (_imod == 1)
			res = b4qi + b4ap1;
		else if (_imod == 2)
			res = b4qi + b4ap2;
		else
			res = b4qi + 0.5 * ( b4ap1 + b4ap2 );

		return res/16.0;
	}


	double P3qg::calcRegular(double x) const
	{
		double nf = static_cast<double>(_nf);
		double YM  = 1.0 / x;
		double Y1  = 1.0 - x;
		double DL  = std::log(x);
		double DL1 = std::log(1.0-x);

		double nf2 = nf * nf;
		double nf3 = nf * nf2;

		// Large-x coefficients
		double x1L5cff = 1.8518519e0 * nf - 4.1152263e-1 * nf2;
		double x1L4cff = 3.5687794e1 * nf - 3.5116598e0 * nf2
			- 8.2304527e-2 * nf3;

		double y1L5cff = 2.8806584e0 * nf + 8.2304527e-1 * nf2;
		double y1L4cff = -4.0511391e1 * nf + 5.5418381e0 * nf2
			+ 1.6460905e-1 * nf3;

		// Small-x coefficients
		double bfkl1   = 3.9357613e3 * nf;

		double x0L6cff = -1.9588477e1 * nf + 2.7654321e0 * nf2;
		double x0L5cff =  2.1573663e1 * nf + 1.7244444e1 * nf2;
		double x0L4cff = -2.8667643e3 * nf + 3.0122403e2 * nf2
			+ 4.1316872e0 * nf3;

		// Base contribution
		double P3QG01 =
			bfkl1 * YM * std::pow(DL, 2)
			+ x0L6cff * std::pow(DL, 6)
			+ x0L5cff * std::pow(DL, 5)
			+ x0L4cff * std::pow(DL, 4)
			+ x1L4cff * std::pow(DL1, 4)
			+ x1L5cff * std::pow(DL1, 5)
			+ y1L4cff * Y1 * std::pow(DL1, 4)
			+ y1L5cff * Y1 * std::pow(DL1, 5);

		double P3qgApp1 = 0.0;
		double P3qgApp2 = 0.0;

		if (nf == 3) {
			P3qgApp1 = P3QG01
				+ 187500.0 * YM * DL
				+ 826060.0 * YM * Y1
				- 150474.0
				+ 226254.0 * x * (2.0 - x)
				+ 577733.0 * DL
				- 180747.0 * std::pow(DL, 2)
				+ 95411.0  * std::pow(DL, 3)
				+ 119.8    * std::pow(DL1, 3)
				+ 7156.3   * std::pow(DL1, 2)
				+ 45790.0  * DL1
				- 95682.0  * DL * DL1;

			P3qgApp2 = P3QG01
				+ 135000.0 * YM * DL
				+ 484742.0 * YM * Y1
				- 11627.0
				- 187478.0 * x * (2.0 - x)
				+ 413512.0 * DL
				- 82500.0  * std::pow(DL, 2)
				+ 29987.0  * std::pow(DL, 3)
				- 850.1    * std::pow(DL1, 3)
				- 11425.0  * std::pow(DL1, 2)
				- 75323.0  * DL1
				+ 282836.0 * DL * DL1;
		}
		else if (nf == 4) {
			P3qgApp1 = P3QG01
				+ 250000.0 * YM * DL
				+ 1089180.0 * YM * Y1
				- 241088.0
				+ 342902.0 * x * (2.0 - x)
				+ 720081.0 * DL
				- 247071.0 * std::pow(DL, 2)
				+ 126405.0 * std::pow(DL, 3)
				+ 272.4    * std::pow(DL1, 3)
				+ 10911.0  * std::pow(DL1, 2)
				+ 60563.0  * DL1
				- 161448.0 * DL * DL1;

			P3qgApp2 = P3QG01
				+ 180000.0 * YM * DL
				+ 634090.0 * YM * Y1
				- 55958.0
				- 208744.0 * x * (2.0 - x)
				+ 501120.0 * DL
				- 116073.0 * std::pow(DL, 2)
				+ 39173.0  * std::pow(DL, 3)
				- 1020.8   * std::pow(DL1, 3)
				- 13864.0  * std::pow(DL1, 2)
				- 100922.0 * DL1
				+ 343243.0 * DL * DL1;
		}
		else if (nf == 5) {
			P3qgApp1 = P3QG01
				+ 312500.0 * YM * DL
				+ 1345700.0 * YM * Y1
				- 350466.0
				+ 480028.0 * x * (2.0 - x)
				+ 837903.0 * DL
				- 315928.0 * std::pow(DL, 2)
				+ 157086.0 * std::pow(DL, 3)
				+ 472.7    * std::pow(DL1, 3)
				+ 15415.0  * std::pow(DL1, 2)
				+ 75644.0  * DL1
				- 244869.0 * DL * DL1;

			P3qgApp2 = P3QG01
				+ 225000.0 * YM * DL
				+ 776837.0 * YM * Y1
				- 119054.0
				- 209530.0 * x * (2.0 - x)
				+ 564202.0 * DL
				- 152181.0 * std::pow(DL, 2)
				+ 48046.0  * std::pow(DL, 3)
				- 1143.8   * std::pow(DL1, 3)
				- 15553.0  * std::pow(DL1, 2)
				- 126212.0 * DL1
				+ 385995.0 * DL * DL1;
		}
		else if (nf == 6) {
			P3qgApp1 = P3QG01
				+ 375000.0 * YM * DL
				+ 1595330.0 * YM * Y1
				- 477729.0
				+ 637552.0 * x * (2.0 - x)
				+ 931556.0 * DL
				- 387017.0 * std::pow(DL, 2)
				+ 187509.0 * std::pow(DL, 3)
				+ 715.5    * std::pow(DL1, 3)
				+ 20710.0  * std::pow(DL1, 2)
				+ 91373.0  * DL1
				- 346374.0 * DL * DL1;

			P3qgApp2 = P3QG01
				+ 270000.0 * YM * DL
				+ 912695.0 * YM * Y1
				- 200034.0
				- 189918.0 * x * (2.0 - x)
				+ 603114.0 * DL
				- 190521.0 * std::pow(DL, 2)
				+ 56661.0  * std::pow(DL, 3)
				- 1224.3   * std::pow(DL1, 3)
				- 16453.0  * std::pow(DL1, 2)
				- 150856.0 * DL1
				+ 410661.0 * DL * DL1;
		}
		else {
			throw std::invalid_argument("Error in P3QGA: invalid nf (must be 3–6)");
		}

		double res{};
		if (_imod == 1)
			res = P3qgApp1;
		else if (_imod == 2)
			res = P3qgApp2;
		else
			res = 0.5 * (P3qgApp1 + P3qgApp2);
		return res/16.0;
	}

	double P3gq::calcRegular(double x) const
	{
		double nf = static_cast<double>(_nf);
		
		double YM  = 1.0 / x;
		double Y1  = 1.0 - x;
		double DL  = std::log(x);
		double DL1 = std::log(Y1);

		int nf2 = nf * nf;
		int nf3 = nf * nf2;   // (kept for strict faithfulness)

		// --------------------------------------------------------------
		// Known large-x coefficients
		// --------------------------------------------------------------

		double x1L5cff = 1.3443073e1 - 5.4869684e-1 * nf;
		double x1L4cff = 3.7539831e2 - 3.4494742e1 * nf
			+ 8.7791495e-1 * nf2;

		double y1L5cff = 2.2222222e1 - 5.4869684e-1 * nf;
		double y1L4cff = 6.6242163e2 - 4.7992684e1 * nf
			+ 8.7791495e-1 * nf2;

		// --------------------------------------------------------------
		// Small-x x^-1 coefficients
		// --------------------------------------------------------------

		double bfkl0 = -8.3086173e3 / 2.25;
		double bfkl1 = (-1.0691199e5 - nf * 9.9638304e2) / 2.25;

		// --------------------------------------------------------------
		// Small-x double logs
		// --------------------------------------------------------------

		double x0L6cff =  5.2235940e1 - 7.3744856e0 * nf;
		double x0L5cff = -2.9221399e2 + 1.8436214e0 * nf;
		double x0L4cff =  7.3106077e3 - 3.7887135e2 * nf
			- 3.2438957e1 * nf2;

		// --------------------------------------------------------------
		// Base function: P3gq01
		// --------------------------------------------------------------

		double P3gq01 =
			bfkl0 * YM * std::pow(DL, 3)
			+ bfkl1 * YM * std::pow(DL, 2)
			+ x0L6cff * std::pow(DL, 6)
			+ x0L5cff * std::pow(DL, 5)
			+ x0L4cff * std::pow(DL, 4)
			+ x1L4cff * std::pow(DL1, 4)
			+ x1L5cff * std::pow(DL1, 5)
			+ y1L4cff * Y1 * std::pow(DL1, 4)
			+ y1L5cff * Y1 * std::pow(DL1, 5);

		double P3gqApp1 = 0.0;
		double P3gqApp2 = 0.0;

		// --------------------------------------------------------------
		// nf-specific approximations
		// --------------------------------------------------------------

		if (_nf == 3) {
			P3gqApp1 =
				P3gq01
				+ 3.5 * bfkl1 * YM * DL
				- 27891.  * YM * Y1
				- 309124.
				+ 1056866. * x * (2.0 - x)
				- 124735. * DL
				- 16246.  * std::pow(DL,2)
				+ 131175. * std::pow(DL,3)
				+ 4970.1  * std::pow(DL1,3)
				+ 60041.  * std::pow(DL1,2)
				+ 343181. * DL1
				- 958330. * DL * DL1;

			P3gqApp2 =
				P3gq01
				+ 7.0 * bfkl1 * YM * DL
				- 1139334. * YM * Y1
				+ 143008.
				- 290390. * x * (2.0 - x)
				- 659492. * DL
				+ 303685. * std::pow(DL,2)
				- 81867.  * std::pow(DL,3)
				+ 1811.8  * std::pow(DL1,3)
				- 465.9   * std::pow(DL1,2)
				- 51206.  * DL1
				+ 274249. * DL * DL1;
		} else if (_nf == 4) {
			P3gqApp1 =
				P3gq01
				+ 3.5 * bfkl1 * YM * DL
				- 8302.8 * YM * Y1
				- 347706.
				+ 1105306. * x * (2.0 - x)
				- 127650. * DL
				- 29728.  * std::pow(DL,2)
				+ 137537. * std::pow(DL,3)
				+ 4658.1  * std::pow(DL1,3)
				+ 59205.  * std::pow(DL1,2)
				+ 345513. * DL1
				- 995120. * DL * DL1;

			P3gqApp2 =
				P3gq01
				+ 7.0 * bfkl1 * YM * DL
				- 1129822. * YM * Y1
				+ 108527.
				- 254166. * x * (2.0 - x)
				- 667254. * DL
				+ 293099. * std::pow(DL,2)
				- 77437.  * std::pow(DL,3)
				+ 1471.3  * std::pow(DL1,3)
				- 1850.3  * std::pow(DL1,2)
				- 52451.  * DL1
				+ 248634. * DL * DL1;
		} else if (_nf == 5) {
			P3gqApp1 =
				P3gq01
				+ 3.5 * bfkl1 * YM * DL
				+ 14035. * YM * Y1
				- 384003.
				+ 1152711. * x * (2.0 - x)
				- 126346. * DL
				- 42967.  * std::pow(DL,2)
				+ 144270. * std::pow(DL,3)
				+ 4385.5  * std::pow(DL1,3)
				+ 58688.  * std::pow(DL1,2)
				+ 348988. * DL1
				- 1031165. * DL * DL1;

			P3gqApp2 =
				P3gq01
				+ 7.0 * bfkl1 * YM * DL
				- 1117561. * YM * Y1
				+ 76329.
				- 218973. * x * (2.0 - x)
				- 670799. * DL
				+ 282763. * std::pow(DL,2)
				- 72633.  * std::pow(DL,3)
				+ 1170.0  * std::pow(DL1,3)
				- 2915.5  * std::pow(DL1,2)
				- 52548.  * DL1
				+ 223771. * DL * DL1;
		} else if (_nf == 6) {
			P3gqApp1 =
				P3gq01
				+ 3.5 * bfkl1 * YM * DL
				+ 39203. * YM * Y1
				- 417914.
				+ 1199042. * x * (2.0 - x)
				- 120750. * DL
				- 55941.  * std::pow(DL,2)
				+ 151383. * std::pow(DL,3)
				+ 4149.2  * std::pow(DL1,3)
				+ 58466.  * std::pow(DL1,2)
				+ 353589. * DL1
				- 1066510. * DL * DL1;

			P3gqApp2 =
				P3gq01
				+ 7.0 * bfkl1 * YM * DL
				- 1102470. * YM * Y1
				+ 46517.
				- 184858. * x * (2.0 - x)
				- 670056. * DL
				+ 272689. * std::pow(DL,2)
				- 67453.  * std::pow(DL,3)
				+ 905.0   * std::pow(DL1,3)
				- 3686.2  * std::pow(DL1,2)
				- 51523.  * DL1
				+ 199594. * DL * DL1;
		} else {
			log(LOG_ERROR, "P3gq::calcRegular()", "nf ({}) out of bounds, not between 3 and 6.", _nf);
		}

		// --------------------------------------------------------------
		// return according to IMOD
		// --------------------------------------------------------------

		double res{};
		if (_imod == 1)
			res = P3gqApp1;
		else if (_imod == 2)
			res = P3gqApp2;
		else
			res = 0.5 * (P3gqApp1 + P3gqApp2);
		return res/16.0;
	}

	double P3gg::calcRegular(double x) const
	{
		double nf = static_cast<double>(_nf);
		
		double YM  = 1.0 / x;
		double Y1  = 1.0 - x;
		double DL  = std::log(x);
		double DL1 = std::log1p(-x);

		double nf2 = nf * nf;
		double nf3 = nf * nf2;

		// Large-x coefficients
		double A4gluon =  40880.330
			- 11714.246  * nf
			+   440.04876 * nf2
			+     7.3627750 * nf3;

		double Ccoeff = 8.5814120e4 - 1.3880515e4 * nf + 1.3511111e2 * nf2;
		double Dcoeff = 5.4482808e4 - 4.3411337e3 * nf - 2.1333333e1 * nf2;

		// 1/x * ln^a(x) terms
		double bfkl0 = -8.308617314e3;
		double bfkl1 = -1.069119905e5 - 9.963830436e2 * nf;

		// Base part
		double P3gg01 =
			bfkl0 * std::pow(DL, 3) * YM
			+ bfkl1 * std::pow(DL, 2) * YM
			+ Ccoeff * DL1
			+ Dcoeff - A4gluon;

		double P3ggApp1 = 0.0;
		double P3ggApp2 = 0.0;

		if (nf == 3) {
			P3ggApp1 = P3gg01
				+ 3.4 * bfkl1 * DL * YM
				- 345063.0 * Y1 * YM
				+ 86650.0 * (1.0 + x * x) * Y1
				+ 158160.0 * DL
				- 15741.0 * Y1 * std::pow(DL1, 2)
				- 9417.0 * Y1 * std::pow(DL1, 3);

			P3ggApp2 = P3gg01
				+ 5.4 * bfkl1 * DL * YM
				- 1265632.0 * Y1 * YM
				- 656644.0 * (1.0 + x * x) * Y1
				- 1352233.0 * DL
				+ 203298.0 * Y1 * std::pow(DL1, 2)
				+ 39112.0 * Y1 * std::pow(DL1, 3);
		}
		else if (nf == 4) {
			P3ggApp1 = P3gg01
				+ 3.4 * bfkl1 * DL * YM
				- 342625.0 * Y1 * YM
				+ 100372.0 * (1.0 + x * x) * Y1
				+ 189167.0 * DL
				- 29762.0 * Y1 * std::pow(DL1, 2)
				- 12102.0 * Y1 * std::pow(DL1, 3);

			P3ggApp2 = P3gg01
				+ 5.4 * bfkl1 * DL * YM
				- 1271540.0 * Y1 * YM
				- 649661.0 * (1.0 + x * x) * Y1
				- 1334919.0 * DL
				+ 191263.0 * Y1 * std::pow(DL1, 2)
				+ 36867.0 * Y1 * std::pow(DL1, 3);
		}
		else if (nf == 5) {
			P3ggApp1 = P3gg01
				+ 3.4 * bfkl1 * DL * YM
				- 337540.0 * Y1 * YM
				+ 119366.0 * (1.0 + x * x) * Y1
				+ 223769.0 * DL
				- 45129.0 * Y1 * std::pow(DL1, 2)
				- 15046.0 * Y1 * std::pow(DL1, 3);

			P3ggApp2 = P3gg01
				+ 5.4 * bfkl1 * DL * YM
				- 1274800.0 * Y1 * YM
				- 637406.0 * (1.0 + x * x) * Y1
				- 1314010.0 * DL
				+ 177882.0 * Y1 * std::pow(DL1, 2)
				+ 34362.0 * Y1 * std::pow(DL1, 3);
		}
		else {
			throw std::invalid_argument("Error in P3GGA: invalid nf (must be 3,4,5)");
		}

		double res{};
		if (_imod == 1)
			res = P3ggApp1;
		else if (_imod == 2)
			res = P3ggApp2;
		else
			res = 0.5 * (P3ggApp1 + P3ggApp2);
		return res/16.0;
	}

	double P3gg::calcPlus(double x) const
	{
		UNUSED(x);
		const double Nf = static_cast<double>(_nf);
		
		const double res = 40880.330 - 11714.246 * Nf + 440.04876 * std::pow(Nf, 2) + 7.3627750 * std::pow(Nf, 3);
		return res/16.0;
	}
	double P3gg::calcDelta(double x) const
	{
		UNUSED(x);
		const double Nf = static_cast<double>(_nf);
		
		const double res = 68587.64 - 18143.983e0 * Nf + 423.81135e0 * std::pow(Nf, 2) + 9.0672154e-1 * std::pow(Nf, 3);
		return res/16.0;
	}

	
	

} // namespace Candia2
