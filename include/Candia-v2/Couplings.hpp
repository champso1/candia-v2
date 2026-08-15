/**
 *  @file Couplings.hpp
 *  @brief Contains the @a AlphaS and @a AlphaQED classes
 */

#pragma once

#include "Candia-v2/Common.hpp"

#include <array>
#include <cmath>

namespace Candia2
{
	class Coupling
	{
	protected:
		uint _order{}; //!< Perturbative order.
		double _Q0{}, _Qf{}; //!< initial and final evolution energies
		double _alpha0{}; //!< Initial value of couplin gat Q0

	    double _beta0{}; //!< \f$\beta_0\f$
		double _beta1{}; //!< \f$\beta_1\f$
		double _beta2{}; //!< \f$\beta_2\f$
		double _beta3{}; //!< \f$\beta_3\f$

		double _mur2_muf2; //!< \f$\frac{mu_R^2}/{mu_F^2}\f$, NO LOG
		double _L{}; //!< \f$\log(\frac{mu_R^2}/{mu_F^2})\f$

		uint _nf{}, _nl{};
		
	public:
		Coupling(uint order, double Q0, double Qf, double alpha0, double mur2_muf2)
			: _order{order}, _Q0{Q0}, _Qf{Qf}, _alpha0{alpha0}, _mur2_muf2{mur2_muf2}, _L{std::log(mur2_muf2)}
		{}
	protected:
		virtual ~Coupling() = default;
		
	public:
		inline virtual double beta0() const { return _beta0; }
		inline virtual double beta1() const { return _beta1; }
		inline virtual double beta2() const { return _beta2; }
		inline virtual double beta3() const { return _beta3; }

		inline virtual double Q0() const { return _Q0; }
		inline virtual double Qf() const { return _Qf; }

		/** @brief evaluates the coupling at Qf given Q0
		 *
		 *  @param Q0 the initial value at which the coupling is evaluated at
		 *  @param Qf the final value to evaluate the coupling at
		 *  @param alpha0 the value of the coupling at the initial energy @a Q0
		 */
		virtual double evaluate(double Q0, double Qf, double alpha0) const = 0;

		/**
		 *  this method will just return a0 and evalutate(q0, qf, a0)
		 */
		virtual std::pair<double,double> initFinalAlpha() const
		{
			return std::make_pair(_alpha0, evaluate(Q0(), Qf(), _alpha0));
		}

	protected:
		virtual double calcBeta0(uint nf, uint nl) const = 0; //!< calculates \f$\beta_0\f$ for @a nf flavors and @a nf leptons
		virtual double calcBeta1(uint nf, uint nl) const = 0; //!< calculates \f$\beta_1\f$ for @a nf flavors and @a nf leptons
		virtual double calcBeta2(uint nf, uint nl) const = 0; //!< calculates \f$\beta_2\f$ for @a nf flavors and @a nf leptons
		virtual double calcBeta3(uint nf, uint nl) const = 0; //!< calculates \f$\beta_3\f$ for @a nf flavors and @a nf leptons

	public:
		/**
		 *  @brief Updates the values of the beta-coefficients at @a nf
		 *  @param nf requested value of \f$n_f\f$
		 */
		virtual void update(uint nf, uint nl)
		{
			_nf = nf;
			_nl = nf;
			_beta0 = calcBeta0(nf, nl);
			_beta1 = calcBeta1(nf, nl);
			_beta2 = calcBeta2(nf, nl);
			_beta3 = calcBeta3(nf, nl);
		}
	};



	class AlphaQED final : public Coupling
	{
	public:
		using Coupling::Coupling;

		inline double evaluate(double Q0, double Qf, double alpha0) const override
		{
			double log_arg = (Qf/Q0)*(Qf/Q0);
		    return alpha0/(1.0 + (alpha0/PI_4)*beta0()*std::log(log_arg));
		}

	private:
		inline double calcBeta0(uint nf, uint nl) const override
		{
			auto nn = nf + nl;
			double fac = nl;
			for (uint i=0; i<nf; ++i)
				fac += NC*Q_QUARK[i]*Q_QUARK[i];
			return -4.0/3.0*fac;
		}
		inline double calcBeta1(uint nf, uint nl) const override
		{
			throw std::runtime_error("TODO: beta1/2/3 for alphaqed");
		}
		inline double calcBeta2(uint nf, uint nl) const override
		{
			throw std::runtime_error("TODO: beta1/2/3 for alphaqed");
		}
		inline double calcBeta3(uint nf, uint nl) const override
		{
			throw std::runtime_error("TODO: beta1/2/3 for alphaqed");
		}

	public:
		void update(uint nf, uint nl) override
		{
			_nf = nf;
			_nl = nf;
			_beta0 = calcBeta0(nf, nl);
			// _beta1 = calcBeta1(nf, nl);
			// _beta2 = calcBeta2(nf, nl);
			// _beta3 = calcBeta3(nf, nl);
		}
	};
}
