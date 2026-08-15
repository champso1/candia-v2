#pragma once

#include "Candia-v2/Common.hpp"
#include "Candia-v2/SplittingFn.hpp"

namespace Candia2
{
	class SplitFuncQED : public SplittingFunction
	{
	public:
		enum class FermionType
		{
			QUARK = 0,
			LEPTON = 1,
			NUM_FERMION_TYPES
		};
	protected:
		using nc_type = std::array<double, static_cast<uint>(FermionType::NUM_FERMION_TYPES)>;
		static nc_type _nc;
		static uint _nl;
		
		FermionType _fermion_type;
		double _qf2{};
	public:
		using SplittingFunction::SplittingFunction;

		inline virtual void set_fermion(FermionType type, double charge)
		{
			_fermion_type = type;
			_qf2 = charge*charge;
		}

		inline static void update(uint nf, double beta0, double log_muf2_mur2, uint nl)
		{
			SplittingFunction::update(nf, beta0, log_muf2_mur2);
			_nl = nl;
		}

		inline virtual double calcPlus(double x) const
		{
			return 0.0;
		}
	};

	class P0ff final: public SplitFuncQED
	{
	public:
		using SplitFuncQED::SplitFuncQED;
		
		inline double calcPlus(double x) const override
		{
			return _qf2*(1.0+x*x);
		}
		inline double calcDelta() const override
		{
			return _qf2*(3.0/2.0);
		}
	};

	class P0fy final: public SplitFuncQED
	{
	public:
		using SplitFuncQED::SplitFuncQED;
		
		inline double calcRegular(double x) const override
		{
			double x2 = 1.0-x*x;
			return _qf2*(1.0 + x2*x2)/x;
		}
	};

	class P0yf final: public SplitFuncQED
	{
	public:
		using SplitFuncQED::SplitFuncQED;
		
		inline double calcRegular(double x) const override
		{
		    double x2 = 1.0-x*x;
			auto nc = _nc[static_cast<uint>(_fermion_type)];
			return nc*_qf2*(x*x + x2*x2);
		}
	};

	class P0yy final: public SplitFuncQED
	{
	public:
		using SplitFuncQED::SplitFuncQED;
		
		inline double calcDelta() const override
		{
			double fac = _nl;
			for (uint nf=0; nf<_nf; ++nf)
				fac +=3*(Q_QUARK[nf]*Q_QUARK[nf]);
			return -(2.0/3.0)*fac;
		}
	};
	
}
