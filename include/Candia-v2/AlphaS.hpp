#ifndef __ALPHAS_HPP
#define __ALPHAS_HPP

#include <array>

#include "Candia-v2/Common.hpp"

namespace Candia2
{

	class AlphaS
	{
	private:
		uint _order{}; //!< Perturbative order.
		double _Q0{}, _alpha0{}; //!< Initial value of \f$\alpha_s\f$ at Q0

	    double _beta0{}, _beta1{}, _beta2{}, _beta3{};
		double _L{}; //!< log (mu_R/mu_F)

		std::array<double, 8> _masses{}; //!< Values of quark masses.
		std::array<double, 8> _pre{}, _post{};
	public:
		AlphaS() = default;
	    AlphaS(
			uint order, double Q0, double alpha0,
			std::array<double,8> const& masses,
			double log_mur2_muf2)
		  : _order{order}, _Q0{Q0}, _alpha0{alpha0}, _L{log_mur2_muf2}, _masses{masses}
		{}

		inline uint order() const { return _order; }
		double masses(uint nf) const;
		uint nff(uint nfi, double Qf);
	    
		double betaFn(double alpha) const;

		inline double beta0() const { return _beta0; };
		inline double beta1() const { return _beta1; };
		inline double beta2() const { return _beta2; };
		inline double beta3() const { return _beta3; };
	public:

		void calculateThresholdValues(double Qf);
		double pre(uint nf) const;
		double post(uint nf) const;

		double evaluate(
			double Q0, double Qf,
			double alpha0) const;

		void update(uint nf);

	private:
		double postMatch(double alpha, uint nf);
		double preMatch(double alpha, uint nf);

		double calcBeta0(uint nf) const;
		double calcBeta1(uint nf) const;
		double calcBeta2(uint nf) const;
		double calcBeta3(uint nf) const;
	};

} // class AlphaS

#endif // __ALPHAS_HPP
