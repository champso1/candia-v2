#include "Candia-v2/SplittingFnQED.hpp"

#include "Candia-v2/SpecialFuncs.hpp"

namespace Candia2
{
	uint SplitFuncQED::_nl;
	double SplitFuncQED::_totalchargefac;

	double SplitFuncQED::pqq(double x) const
	{
		return 1+x*x; // coefficient of the plus distribution
	}
	double SplitFuncQED::pqg(double x) const
	{
		double x1 = 1.0-x;
		return x*x + x1*x1;
	}
	double SplitFuncQED::pgq(double x) const
	{
		double x1 = 1.0-x;
		return (1.0 + x1*x1)/x;
	}
	double SplitFuncQED::ps(double x) const
	{
		double lx = std::log(x);
		return 20.0/(9.0*x) - 2.0 + 6.0*x - (56.0/9.0)*x*x + (1 + 5.0*x + (8.0/3.0)*x*x)*lx - (1.0+x)*lx*lx;
	}
	double SplitFuncQED::S2(double x) const
	{
		double lx = std::log(x);
		double lx2 = lx*lx;
		double lpx = std::log1p(x);

		return lx2/2.0 - Zeta2 - 2*Li2(-x) - 2.0*lx*lpx;
	}


	double P1ffV::regular_nofac1(double x) const
	{
		double x1 = 1.0-x;
		double lx = std::log(x);
		double lx2 = lx*lx;
		double lx1 = std::log(x1);

		return ((3.0 + 7.0*x)/2.0)*lx + ((1.0+x)/2.0)*lx2 + 5.0*x1;
	}
	double P1ffV::regular_nofac2(double x) const
	{
		return (4.0/3.0)*(1.0-x);
	}
	double P1ffV::plus_nofac1(double x) const
	{
		double x1 = 1.0-x;
		double lx = std::log(x);
		double lx1 = std::log(x1);
		return (2.0*lx*lx1 + (3.0/2.0)*lx)*pqq(x);
	}
	double P1ffV::plus_nofac2(double x) const
	{
		double lx = std::log(x);
		return ((2.0/3.0)*lx + 10.0/9.0)*pqq(x);
	}
	double P1ffV::delta_nofac1() const
	{
		return PI_2/2.0 - 3.0/8.0 - 6.0*Zeta3;
	}
	double P1ffV::delta_nofac2() const
	{
		return 2.0*PI_2/9.0 + 1.0/6.0;
	}



	double P1fy::regular_nofac(double x) const
	{
		double dx = (1.0-x)/x;
		double lx = std::log(x);
		double lx2 = lx*lx;
		double lx1 = std::log1p(-x);
		double ldx = std::log(dx);
		double ldx2 = ldx*ldx;
			
		return 0.5*(4.0 - 9.0*x - (1.0 - 4.0*x)*lx - (1.0 - 2.0*x)*lx2 + 4.0*lx1 + pqg(x)*(2.0*ldx2 - 4.0*ldx - 2.0*PI_2/3.0 + 10));
	}



	double P1yf::regular_nofac1(double x) const
	{
		double lx = std::log(x);
		double lx2 = lx*lx;
		double lx1 = std::log1p(-x);
		double lx12 = lx1*lx1;
		return -(3.0*lx1 + lx12)*pgq(x) + (2.0 + (7.0/2.0)*x)*lx - (1.0 - x/2.0)*lx2 - 2.0*x*lx1 - (7.0/2.0)*x - 5.0/2.0;
	}
	double P1yf::regular_nofac2(double x) const
	{
		double lx1 = std::log1p(-x);
		return (4.0/3.0)*x + pgq(x)*(20.0/9.0 + (4.0/3.0)*lx1);
	}



	double P11qB::regular_nofac(double x) const
	{
		double lx = std::log(x);
		double lx2 = lx*lx;
		double lx1 = std::log1p(-x);
		double dx = (1.0-x)/x;
		double ldx = std::log(dx);
		double ldx2 = ldx*ldx;
		return 0.5*(4.0 - 9.0*x - (1.0 - 4.0*x)*lx - (1.0 - 2.0*x)*lx2 + 4.0*lx1 + pqg(x)*(2.0*ldx2 - 4.0*ldx - 2*PI_2/3.0 + 10.0));
	}



	double P11bB::regular_nofac(double x) const
	{
		double x2 = x*x;
		double lx = std::log(x);
		double lx2 = lx*lx;
		return -16.0 + 8.0*x + (20.0/3.0)*x2 + 4.0/(3.0*x) - (6.0 + 10.0*x)*lx - 2.0*(1.0+x)*lx2;
	}




	double P11qqV::regular_nofac(double x) const
	{
		double lx = std::log(x);
		double lx2 = lx*lx;
		return -2.0*(((3.0 + 7.0*x)/2)*lx + ((1.0 + x)/2.0)*lx2 + 5.0*(1.0-x));
	}
	double P11qqV::plus_nofac(double x) const
	{
		double lx = std::log(x);
		double lx1 = std::log1p(-x);
		return -2.0*(2.0*lx1 + 3.0/2.0)*lx*pqq(x);
	}
	double P11qqV::delta_nofac() const
	{
		return -2.0*(PI_2/2.0 - 3.0/8.0 - 6.0*Zeta3);
	}



	double P11qqbarV::regular_nofac(double x) const
	{
		double lx = std::log(x);
		return 2.0*(4.0*(1.0-x) + 2*(1.0+x)*lx);
	}
	double P11qqbarV::plus_nofac(double x) const
	{
		double lx = std::log(x);
		double lx1 = std::log1p(-x);
		return 2.0*(2.0*pqq(-x)*S2(x));
	}



	double P11bq::regular_nofac(double x) const
	{
		double lx1 = std::log1p(-x);
		double lx12 = lx1*lx1;
		double lx = std::log(x);
		double lx2 = lx*lx;
		return -(3.0*lx1 + lx12)*pgq(x) + (2.0 + (7.0/2.0)*x)*lx - (1.0 - x/2.0)*lx2 - 2.0*x*lx1 - (7.0/2.0)*x - 5.0/2.0;
	}




	
}
