#include "Candia-v2/Common.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"
#include "Candia-v2/SpecialFuncs.hpp"

namespace Candia2
{

	double A2ns::Regular(const double x) const
	{
		double L = std::log(x);

		return CF*TR*((1.0 + x*x)/(1.0-x)*((2.0/3.0)*L*L + (20.0/9.0)*L) 
					  + (8.0/3.0)*(1.0-x)*L + 44.0/27.0 - (268.0/27.0)*x);
	}

	double A2ns::Plus(const double x) const
	{
		UNUSED(x);
		return CF*TR*(224.0/27.0);
	}

	double A2ns::Delta(const double x) const
	{
		UNUSED(x);
		return CF*TR*((-8.0/3.0)*Zeta3 + (40.0/9.0)*Zeta2 + 73.0/18.0);
	}


	double A2gq::Regular(const double x) const
	{
		double M = std::log1p(-x);
		return CF*TR*((4.0/3.0)*(2.0/x - 2.0 + x)*M*M
					  + (8.0/9.0)*(10.0/x - 10.0 + 8.0*x)*M
					  + (448.0/x - 448.0 + 344.0*x)/27.0);
	}

	double A2gq::Plus(const double x) const
	{
		UNUSED(x);
		return 0.0;
	}

	double A2gq::Delta(const double x) const
	{
		UNUSED(x);
		return 0.0;
	}



    double A2gg::Regular(const double x) const
	{
		double L = std::log(x);
		double M = std::log1p(-x);

		return CF*TR*((4.0/3.0)*(1.0+x)*L*L*L
					  + (6.0 + 10.0*x)*L*L 
					  + (32.0 + 48.0*x)*L 
					  - 8.0/x + 80.0 - 48.0*x - 24.0*x*x) 
				+ NC*TR*(4.0/3.0*(1.0+x)*L*L 
				+ (52.0 + 88.0*x)*L/9.0 - 4.0/3.0*x*M
				+ (556.0/x - 628.0 + 548.0*x - 700.0*x*x)/27.0);
	}

	double A2gg::Plus(const double x) const
	{
		UNUSED(x);
		return NC*TR*224.0/27.0;
	}

	double A2gg::Delta(const double x) const
	{
		UNUSED(x);
		return -CF*TR*15.0 + NC*TR*10.0/(9.0*27.0);
	}




	double A2hq::Regular(const double x) const
	{
		double L = std::log(x);


		return CF*TR*((1.0+x)*(32.0*S12(1.0-x) + 16.0*L*Li2(1.0-x) - 16.0*Zeta2*L - 4.0/3.0*L*L*L) + (32.0/(3.0*x) + 8.0 - 8.0*x - 32.0/3.0*x*x)*Li2(1.0-x) + (-32.0/(3.0*x) - 8.0 + 8.0*x + 32.0/3.0*x*x)*Zeta2 + (2.0 + 10.0*x + 16.0/3.0*x*x)*L*L - (56.0/3.0 + 88.0/3.0*x + 448.0/9.0*x*x)*L - 448.0/(27.0*x) - 4.0/3.0 - 124.0/3.0*x + 1600.0/27.0*x*x);
	}

	double A2hq::Plus(const double x) const
	{
		UNUSED(x);
		return 0.0;
	}

	double A2hq::Delta(const double x) const
	{
		UNUSED(x);
		return 0.0;
	}




    double A2hg::Regular(const double x) const
	{
		double L = std::log(x);
		double M = std::log1p(-x);
		double P = std::log1p(x);
		double S = S12(1.0-x);

		return CF*TR*((1.0 - 2.0*x + 2.0*x*x)*(8.0*Zeta3 + 4.0/3.0*M*M*M - 8.0*M*Li2(1.0-x) + 8.0*Zeta2*L - 4.0*L*M*M + 2.0/3.0*L*L*L - 8.0*L*Li2(1.0-x) + 8.0*Li3(1.0-x) - 24.0*S) + x*x*(-16.0*Zeta2*L + 4.0/3.0*L*L*L + 16.0*L*Li2(1.0-x) + 32.0*S) - (4.0 + 96.0*x - 64.0*x*x)*Li2(1.0-x) - (4.0 - 48.0*x + 40.0*x*x)*Zeta2 - (8.0 + 48.0*x - 24.0*x*x)*L*M + (4.0 + 8.0*x - 12.0*x*x)*M*M - (1.0 + 12.0*x - 20.0*x*x)*L*L - (52.0*x - 48.0*x*x)*M - (16.0 + 18.0*x + 48.0*x*x)*L + 26.0 - 82.0*x + 80.0*x*x) + NC*TR*((1.0 - 2.0*x + 2.0*x*x)*(-4.0/3.0*M*M*M + 8.0*M*Li2(1.0-x) - 8.0*Li3(1.0-x))+(1.0 + 2.0*x + 2.0*x*x)*(-8.0*Zeta2*P - 16.0*P*Li2(-x) - 8.0*L*P*P + 4.0*L*L*P + 8.0*L*Li2(-x) - 8.0*Li3(-x) - 16.0*S12(-x)) + (16.0 + 64.0*x)*(2.0*S + L*Li2(1.0-x)) - (4.0/3.0 + 8.0/3.0*x)*L*L*L + (8.0 - 32.0*x + 16.0*x*x)*Zeta3 - (16.0 + 64.0*x)*Zeta2*L + (16.0*x + 16.0*x*x)*(Li2(-x) + L*P) + (32.0/(3.0*x) + 12.0 + 64.0*x - 272.0/3.0*x*x)*Li2(1.0-x) - (12.0 + 48.0*x - 260.0/3.0*x*x + 32.0/(3.0*x))*Zeta2 - 4.0*x*x*L*M - (2.0 + 8.0*x - 10.0*x*x)*M*M + (2.0 + 8.0*x + 46.0/3.0*x*x)*L*L + (4.0 + 16.0*x - 16.0*x*x)*M - (56.0/3.0 + 172.0/3.0*x + 1600.0/9.0*x*x)*L - 448.0/(27.0*x) - 4.0/3.0 - 628.0/3.0*x + 6352.0/27.0*x*x);
	}

	double A2hg::Plus(const double x) const
	{
		UNUSED(x);
		return 0.0;
	}

	double A2hg::Delta(const double x) const
	{
		UNUSED(x);
		return 0.0;
	}
	
};
