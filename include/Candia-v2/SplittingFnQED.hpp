#pragma once

#include "Candia-v2/Common.hpp"
#include "Candia-v2/SplittingFn.hpp"
#include "Candia-v2/SpecialFuncs.hpp"

#include <cmath>

namespace Candia2
{
    class SplitFuncQED : public SplittingFunction
	{
	protected:
		static uint _nl;
		static double _totalchargefac;
	public:
		using SplittingFunction::SplittingFunction;

		inline static void update(uint nf, double beta0, double log_muf2_mur2, uint nl)
		{
			SplittingFunction::update(nf, beta0, log_muf2_mur2);
			_nl = nl;

			_totalchargefac = _nl;
			for (uint i=0; i<_nf; ++i)
				_totalchargefac += NC*Q_QUARK[i];
		}

	protected:
		inline double pqq(double x) const
		{
			return 1+x*x; // coefficient of the plus distribution
		}
		inline double pqg(double x) const
		{
			double x1 = 1.0-x;
			return x*x + x1*x1;
		}
		inline double pgq(double x) const
		{
			double x1 = 1.0-x;
			return (1.0 + x1*x1)/x;
		}
		inline double ps(double x) const
		{
			double lx = std::log(x);
			return 20.0/(9.0*x) - 2.0 + 6.0*x - (56.0/9.0)*x*x + (1 + 5.0*x + (8.0/3.0)*x*x)*lx - (1.0+x)*lx*lx;
		}
		inline double S2(double x) const
		{
			double lx = std::log(x);
			double lx2 = lx*lx;
			double lpx = std::log1p(x);

			return lx2/2.0 - Zeta2 - 2*Li2(-x) - 2.0*lx*lpx;
		}
	};

	class P0ff : public SplitFuncQED
	{
	private:
		double _fac;
	public:
		P0ff(double fac) : _fac{fac} {}
	private:
		inline double plus_nofac(double x) const
		{
			return pqq(x);
		}
		inline double delta_nofac() const
		{
			return 3.0/2.0;
		}
	public:
		inline double calcPlus(double x) const override { return _fac*plus_nofac(x); }
		inline double calcDelta() const override { return _fac*delta_nofac(); }
	};

	class P0uu final : public P0ff
	{
	public:
		P0uu() : P0ff(Q_QUARK[0]*Q_QUARK[0]) {}
	};
	class P0dd final : public P0ff
	{
	public:
		P0dd() : P0ff(Q_QUARK[1]*Q_QUARK[1]) {}
	};
	class P0ll final : public P0ff
	{
	public:
		P0ll() : P0ff(1.0) {}
	};


	class P0fy : public SplitFuncQED
	{
	private:
		double _fac;
	public:
		P0fy(double fac) : _fac{fac} {}
	private:
		inline double regular_nofac(double x) const
		{
			return pqg(x);
		}
	public:
		inline double calcRegular(double x) const override { return _fac*regular_nofac(x); }
	};

	class P0uy final : public P0fy
	{
	public:
		P0uy() : P0fy(NC*Q_QUARK[0]*Q_QUARK[0]) {}
	};
	class P0dy final : public P0fy
	{
	public:
		P0dy() : P0fy(NC*Q_QUARK[1]*Q_QUARK[1]) {}
	};
	class P0ly final : public P0fy
	{
	public:
		P0ly() : P0fy(1) {}
	};

	

	class P0yf : public SplitFuncQED
	{
	private:
		double _fac;
	public:
		P0yf(double fac) : _fac{fac} {}
	private:
		inline double regular_nofac(double x) const
		{
		    return pgq(x);
		}
	public:
		inline double calcRegular(double x) const override { return _fac*regular_nofac(x); }
	};

	class P0yu final : public P0yf
	{
    public:
		P0yu() : P0yf(Q_QUARK[0]*Q_QUARK[0]) {}
	};
	class P0yd final : public P0yf
	{
    public:
		P0yd() : P0yf(Q_QUARK[1]*Q_QUARK[1]) {}
	};
	class P0yl final : public P0yf
	{
    public:
		P0yl() : P0yf(1.0) {}
	};


	class P0yy final: public SplitFuncQED
	{
	public:
		using SplitFuncQED::SplitFuncQED;
		
		inline double calcDelta() const override
		{
			return _totalchargefac*(-2.0/3.0);
		}
	};
	

	class P1ffV : public SplitFuncQED
	{
	private:
		double _fac1, _fac2;
	public:
		using SplitFuncQED::SplitFuncQED;
		P1ffV(double fac1, double fac2) : _fac1{fac1}, _fac2{fac2} {}
	private:
		inline double regular_nofac1(double x) const
		{
			double x1 = 1.0-x;
			double lx = std::log(x);
			double lx2 = lx*lx;
			double lx1 = std::log(x1);

			return ((3.0 + 7.0*x)/2.0)*lx + ((1.0+x)/2.0)*lx2 + 5.0*x1;
		}
		inline double regular_nofac2(double x) const
		{
			return (4.0/3.0)*(1.0-x);
		}
		inline double plus_nofac1(double x) const
		{
			double x1 = 1.0-x;
			double lx = std::log(x);
			double lx1 = std::log(x1);
			return (2.0*lx*lx1 + (3.0/2.0)*lx)*pqq(x);
		}
		inline double plus_nofac2(double x) const
		{
			double lx = std::log(x);
			return ((2.0/3.0)*lx + 10.0/9.0)*pqq(x);
		}
		inline double delta_nofac1() const
		{
			return PI_2/2.0 - 3.0/8.0 - 6.0*Zeta3;
		}
		inline double delta_nofac2() const
		{
			return 2.0*PI_2/9.0 + 1.0/6.0;
		}
	public:
		inline double calcRegular(double x) const override { return _fac1*regular_nofac1(x) + _totalchargefac*_fac2*regular_nofac2(x); }
		inline double calcPlus(double x) const override { return _fac1*plus_nofac1(x) + _totalchargefac*_fac2*plus_nofac2(x); }
		inline double calcDelta() const override { return _fac1*delta_nofac1() + _totalchargefac*_fac2*delta_nofac2(); }
	};

	class P1uuV final : public P1ffV
	{
	public:
		P1uuV() : P1ffV(-power<4>(Q_QUARK[0]), -power<4>(Q_QUARK[0])) {}
	};
	class P1ddV final : public P1ffV
	{
	public:
		P1ddV() : P1ffV(-power<4>(Q_QUARK[1]), -power<4>(Q_QUARK[1])) {}
	};
	class P1llV final : public P1ffV
	{
	public:
		P1llV() : P1ffV(-1.0, -1.0) {}
	};

	class P1ffbarV : public SplitFuncQED
	{
	private:
		double _fac;
	public:
		P1ffbarV(double fac) : _fac{fac} {}
	private:
		inline double regular_nofac(double x) const
		{
			return 4.0*(1.0-x) + 2.0*(1.0+x)*std::log(x);
		}
		inline double plus_nofac(double x) const
		{
			return 2*pqq(-x)*S2(x);
		}
	public:
		inline double calcRegular(double x) const override { return _fac*regular_nofac(x); }
		inline double calcPlus(double x) const override { return _fac*plus_nofac(x); }
	};

	class P1uubarV final : public P1ffbarV
	{
	public:
		P1uubarV() : P1ffbarV(power<4>(Q_QUARK[0])) {}
	};
	class P1ddbarV final : public P1ffbarV
	{
	public:
		P1ddbarV() : P1ffbarV(power<4>(Q_QUARK[1])) {}
	};
	class P1llbarV final : public P1ffbarV
	{
	public:
		P1llbarV() : P1ffbarV(1) {}
	};

	class P1fFS : public SplitFuncQED
	{
	private:
		double _fac;
	public:
		P1fFS(double fac) : _fac{fac} {}
	private:
		inline double regular_nofac(double x) const { return ps(x); }
	public:
		inline double calcRegular(double x) const override { return _fac*regular_nofac(x); }
	};
	class P1uUS final : public P1fFS
	{
	public:
		P1uUS() : P1fFS(NC*power<2>(Q_QUARK[0])*power<2>(Q_QUARK[0])) {}
	};
	class P1uDS final : public P1fFS
	{
	public:
		P1uDS() : P1fFS(NC*power<2>(Q_QUARK[0])*power<2>(Q_QUARK[1])) {}
	};
	class P1dDS final : public P1fFS
	{
	public:
		P1dDS() : P1fFS(NC*power<2>(Q_QUARK[1])*power<2>(Q_QUARK[1])) {}
	};
	class P1lLS final : public P1fFS
	{
	public:
		P1lLS() : P1fFS(1.0) {}
	};


	class P1fy : public SplitFuncQED
	{
	private:
		double _fac;
	public:
		P1fy(double fac) : _fac{fac} {}
	private:
		inline double regular_nofac(double x) const
		{
			double dx = (1.0-x)/x;
			double lx = std::log(x);
			double lx2 = lx*lx;
			double lx1 = std::log1p(-x);
			double ldx = std::log(dx);
			double ldx2 = ldx*ldx;
			
			return 0.5*(4.0 - 9.0*x - (1.0 - 4.0*x)*lx - (1.0 - 2.0*x)*lx2 + 4.0*lx1 + pqg(x)*(2.0*ldx2 - 4.0*ldx - 2.0*PI_2/3.0 + 10));
		}
	public:
	    inline double calcRegular(double x) const override { return _fac*regular_nofac(x); }
	};
	class P1uy final : public P1fy
	{
	public:
		P1uy() : P1fy(NC*power<4>(Q_QUARK[0])) {}
	};
	class P1dy final : public P1fy
	{
	public:
		P1dy() : P1fy(NC*power<4>(Q_QUARK[1])) {}
	};
	class P1ly final : public P1fy
	{
	public:
		P1ly() : P1fy(1.0) {}
	};

	class P1yf : public SplitFuncQED
	{
	private:
		double _fac1, _fac2;
	public:
		P1yf(double fac1, double fac2) : _fac1{fac1}, _fac2{fac2} {}
	private:
		inline double regular_nofac1(double x) const
		{
			double lx = std::log(x);
			double lx2 = lx*lx;
			double lx1 = std::log1p(-x);
			double lx12 = lx1*lx1;
			return -(3.0*lx1 + lx12)*pgq(x) + (2.0 + (7.0/2.0)*x)*lx - (1.0 - x/2.0)*lx2 - 2.0*x*lx1 - (7.0/2.0)*x - 5.0/2.0;
		}
		inline double regular_nofac2(double x) const
		{
			double lx1 = std::log1p(-x);
			return (4.0/3.0)*x + pgq(x)*(20.0/9.0 + (4.0/3.0)*lx1);
		}
	public:
		inline double calcRegular(double x) const override { return _fac1*_totalchargefac*regular_nofac1(x) - _fac2*_totalchargefac*regular_nofac2(x); }
	};
	class P1yu final : public P1yf
	{
	public:
		P1yu() : P1yf(power<4>(Q_QUARK[0]), power<2>(Q_QUARK[0])) {}
	};
	class P1yd final : public P1yf
	{
	public:
		P1yd() : P1yf(power<4>(Q_QUARK[1]), power<2>(Q_QUARK[1])) {}
	};
	class P1yl final : public P1yf
	{
	public:
		P1yl() : P1yf(1.0,1.0) {}
	};

	class P11qB : public SplitFuncQED
	{
	private:
		double _fac;
	public:
		using SplitFuncQED::SplitFuncQED;
		P11qB(double fac) : _fac{fac} {}
	private:
		inline double regular_nofac(double x) const
		{
			double lx = std::log(x);
			double lx2 = lx*lx;
			double lx1 = std::log1p(-x);
			double dx = (1.0-x)/x;
			double ldx = std::log(dx);
			double ldx2 = ldx*ldx;
			return 0.5*(4.0 - 9.0*x - (1.0 - 4.0*x)*lx - (1.0 - 2.0*x)*lx2 + 4.0*lx1 + pqg(x)*(2.0*ldx2 - 4.0*ldx - 2*PI_2/3.0 + 10.0));
		}
	public:
		inline double calcRegular(double x) const override{ return _fac*regular_nofac(x); }
	};
	class P11uy final : public P11qB
	{
	public:
		P11uy() : P11qB(CF*NC*Q_QUARK[0]*Q_QUARK[0]) {}
	};
	class P11dy final : public P11qB
	{
	public:
		P11dy() : P11qB(CF*NC*Q_QUARK[1]*Q_QUARK[1]) {}
	};
	class P11ug final : public P11qB
	{
	public:
		P11ug() : P11qB(TR*Q_QUARK[0]*Q_QUARK[0]) {}
	};
	class P11dg final : public P11qB
	{
	public:
		P11dg() : P11qB(TR*Q_QUARK[1]*Q_QUARK[1]) {}
	};


	class P11bB : public SplitFuncQED
	{
	private:
		double _fac;
		double _chsq_sum;
	public:
		P11bB(double fac) : _fac{fac}
		{
			_chsq_sum = 0;
			for (uint i=0; i<_nf; ++i)
				_chsq_sum += Q_QUARK[i]*Q_QUARK[i];
		}
	private:
		inline double regular_nofac(double x) const
		{
			double x2 = x*x;
			double lx = std::log(x);
			double lx2 = lx*lx;
			return -16.0 + 8.0*x + (20.0/3.0)*x2 + 4.0/(3.0*x) - (6.0 + 10.0*x)*lx - 2.0*(1.0+x)*lx2;
		}
	public:
		inline double calcRegular(double x) const override { return _fac*_chsq_sum*regular_nofac(x); }
	};
	class P11gy final : public P11bB
	{
	public:
		P11gy() : P11bB(CF*NC) {}
	};
	class P11yg final : public P11bB
	{
	public:
		P11yg() : P11bB(TR) {}
	};

	
	class P11bb : public SplitFuncQED
	{
	private:
		double _fac;
		double _chsq_sum;
	public:
		P11bb(double fac) : _fac{fac}
		{
			_chsq_sum = 0;
			for (uint i=0; i<_nf; ++i)
				_chsq_sum += Q_QUARK[i]*Q_QUARK[i];
		}
	private:
		inline double delta_nofac() const { return 1.0; }
	public:
		inline double calcDelta() const override { return -_fac*_chsq_sum*delta_nofac(); }
	};
	class P11yy final : public P11bb
	{
	public:
		P11yy() : P11bb(CF*NC) {}
	};
	class P11gg final : public P11bb
	{
		P11gg() : P11bb(TR) {}
	};

	class P11qqV : public SplitFuncQED
	{
	private:
		double _fac;
	public:
		P11qqV(double fac) : _fac{fac} {}
	private:
		inline double regular_nofac(double x) const
		{
			double lx = std::log(x);
			double lx2 = lx*lx;
			return -2.0*(((3.0 + 7.0*x)/2)*lx + ((1.0 + x)/2.0)*lx2 + 5.0*(1.0-x));
		}
		inline double plus_nofac(double x) const
		{
			double lx = std::log(x);
			double lx1 = std::log1p(-x);
			return -2.0*(2.0*lx1 + 3.0/2.0)*lx*pqq(x);
		}
		inline double delta_nofac() const
		{
			return -2.0*(PI_2/2.0 - 3.0/8.0 - 6.0*Zeta3);
		}
	public:
		inline double calcRegular(double x) const override { return _fac*regular_nofac(x); }
		inline double calcPlus(double x) const override { return _fac*plus_nofac(x); }
		inline double calcDelta() const override { return _fac*delta_nofac(); }
	};
	class P11uuV final : public P11qqV
	{
	public:
		P11uuV() : P11qqV(CF*Q_QUARK[0]*Q_QUARK[0]) {}
	};
	class P11ddV final : public P11qqV
	{
	public:
		P11ddV() : P11qqV(CF*Q_QUARK[1]*Q_QUARK[1]) {}
	};

	class P11qqbarV : public SplitFuncQED
	{
	private:
		double _fac;
	public:
		P11qqbarV(double fac) : _fac{fac} {}
	private:
		inline double regular_nofac(double x) const
		{
			double lx = std::log(x);
			return 2.0*(4.0*(1.0-x) + 2*(1.0+x)*lx);
		}
		inline double plus_nofac(double x) const
		{
			double lx = std::log(x);
			double lx1 = std::log1p(-x);
			return 2.0*(2.0*pqq(-x)*S2(x));
		}
	public:
		inline double calcRegular(double x) const override { return _fac*regular_nofac(x); }
		inline double calcPlus(double x) const override { return _fac*plus_nofac(x); }
	};
	class P11uubarV final : public P11qqbarV
	{
	public:
		P11uubarV() : P11qqbarV(CF*Q_QUARK[0]*Q_QUARK[0]) {}
	};
	class P11ddbarV final : public P11qqbarV
	{
	public:
		P11ddbarV() : P11qqbarV(CF*Q_QUARK[1]*Q_QUARK[1]) {}
	};

	class P11bq : public SplitFuncQED
	{
	private:
		double _fac;
	public:
		P11bq(double fac) : _fac{fac} {}
	private:
		inline double regular_nofac(double x) const
		{
			double lx1 = std::log1p(-x);
			double lx12 = lx1*lx1;
			double lx = std::log(x);
			double lx2 = lx*lx;
			return -(3.0*lx1 + lx12)*pgq(x) + (2.0 + (7.0/2.0)*x)*lx - (1.0 - x/2.0)*lx2 - 2.0*x*lx1 - (7.0/2.0)*x - 5.0/2.0;
		}
	public:
		inline double calcRegular(double x) const override { return _fac*regular_nofac(x); }
	};
	class P11gu final : public P11bq
	{
	public:
		P11gu() : P11bq(CF*Q_QUARK[0]*Q_QUARK[0]) {}
	};
	class P11gd final : public P11bq
	{
	public:
		P11gd() : P11bq(CF*Q_QUARK[1]*Q_QUARK[1]) {}
	};
	class P11yu final : public P11bq
	{
	public:
		P11yu() : P11bq(CF*Q_QUARK[0]*Q_QUARK[0]) {}
	};
	class P11yd final : public P11bq
	{
	public:
		P11yd() : P11bq(CF*Q_QUARK[1]*Q_QUARK[1]) {}
	};
}
