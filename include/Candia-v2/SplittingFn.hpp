#ifndef __SPLITTINGFN_HPP
#define __SPLITTINGFN_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"

namespace Candia2
{	
	class SplittingFunction : public Expression
	{
	protected:
		static uint _nf;      //!< number of active/currently massless flavors
		static double _beta0; //!< beta0 coefficient (for P0gg)
		static double _kr;    //!< log of mu_f/mu_r

		SplittingFunction() = default;
	public:
		virtual ~SplittingFunction() = default;

		inline static uint nf() { return _nf; }
		inline static double beta0() { return _beta0; }
		inline static double kr() { return _kr; }
		inline static void kr(double kr) { _kr = kr; }
		
		inline static void update(uint nf, double beta0)
		{ _nf = nf; _beta0 = beta0; }
	};



	class P0ns final : public SplittingFunction
	{
	public:
		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P0qq final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P0qg final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
	};

	class P0gq final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
	};

	class P0gg final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};



	class P1nsp final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P1nsm final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P1qq final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P1qg final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
	};

	class P1gq final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
	};

	class P1gg final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};




	class P2nsp final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P2nsm final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P2nsv final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P2ps final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
	};

	class P2qq final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P2qg final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
	};

	class P2gq final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
	};

	class P2gg final: public SplittingFunction
	{
	public:
		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};



	
	class P3nsp final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use
		
	public:
		P3nsp(const uint imod=3) : _imod(imod) { }

		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P3nsm final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use
		
	public:
		P3nsm(const uint imod=3) : _imod(imod)  { }

		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P3nsv final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use
		
	public:
		P3nsv(const uint imod=3) : _imod(imod) { }

		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P3ps final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use

	public:
		P3ps(const uint imod=3) : _imod(imod) { }

		double regular(const double x) const override;
	};

	class P3qq final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use

	public:
		P3qq(const uint imod=3) : _imod(imod) { }

		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};

	class P3qg final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use

	public:
		P3qg(const uint imod=3) : _imod(imod) { }

		double regular(const double x) const override;
	};

	class P3gq final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use

	public:
		P3gq(const uint imod=3) : _imod(imod) { }

		double regular(const double x) const override;
	};

	class P3gg final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use

	public:
		P3gg(const uint imod=3) : _imod(imod) { }

		double regular(const double x) const override;
		double plus(const double x) const override;
		double delta(const double x) const override;
	};
} // namespace Candia2


#endif // __SPLITTINGFN_HPP
