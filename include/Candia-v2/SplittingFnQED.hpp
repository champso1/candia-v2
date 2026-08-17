#pragma once

#include "Candia-v2/Common.hpp"
#include "Candia-v2/SplittingFn.hpp"

namespace Candia2
{
	class SplitFuncQED : public SplittingFunction
	{
	protected:
		static uint _nl;
	public:
		using SplittingFunction::SplittingFunction;

		inline static void update(uint nf, double beta0, double log_muf2_mur2, uint nl)
		{
			SplittingFunction::update(nf, beta0, log_muf2_mur2);
			_nl = nl;
		}
	};

	class P0ff final: public SplitFuncQED
	{
	public:
		using SplitFuncQED::SplitFuncQED;
		
		inline double calcPlus(double x) const override
		{
			return 1.0+x*x;
		}
		inline double calcDelta() const override
		{
			return 3.0/2.0;
		}
	};

	class P0fy final: public SplitFuncQED
	{
	public:
		using SplitFuncQED::SplitFuncQED;
		
		inline double calcRegular(double x) const override
		{
			double x1 = 1.0-x;
			return x*x + x1*x1;
		}
	};

	class P0yf final: public SplitFuncQED
	{
	public:
		using SplitFuncQED::SplitFuncQED;
		
		inline double calcRegular(double x) const override
		{
		    double x1 = 1.0-x;
			return (1.0 + x1*x1)/x;
		}
	};

	class P0yy final: public SplitFuncQED
	{
	public:
		using SplitFuncQED::SplitFuncQED;
		
		inline double calcDelta() const override
		{
			return -2.0/3.0;
		}
	};
	
}
