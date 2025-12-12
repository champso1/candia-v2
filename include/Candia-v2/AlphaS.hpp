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
		uint _nf{};    //!< number of active flavors (set via setFFNS())
		uint _nfi{}, _nff{}; //!< initial and final number of flavors
		double _Q0{}, _Qf{}; //!< initial and final evolution energies
		double _alpha0{}; //!< Initial value of \f$\alpha_s\f$ at Q0

	    double _beta0{}, _beta1{}, _beta2{}, _beta3{};
		double _L{}; //!< log (mu_R/mu_F)

		std::array<double, 8> _masses{}; //!< Values of quark masses.
		std::array<double, 8> _pre{}, _post{};

	    enum Scheme : int
		{
			UNSET = -1,
			FIXED = 0,
			VARIABLE = 1
		};
		int _scheme;
	public:
		AlphaS() = default;
	    AlphaS(uint order, double Q0, double Qf, double alpha0, double log_mur2_muf2)
			: _order{order}, _Q0{Q0}, _Qf{Qf}, _alpha0{alpha0}, _L{log_mur2_muf2}, _scheme{UNSET}
		{}

		void setFFNS(uint nf);
		void setVFNS(std::array<double, 8>, uint nfi);

		double masses(uint nf) const;
		inline uint nfi() const { return _nfi; }
		inline uint nff() const { return _nff; }
	    
		double betaFn(double alpha) const;

		inline double beta0() const { return _beta0; };
		inline double beta1() const { return _beta1; };
		inline double beta2() const { return _beta2; };
		inline double beta3() const { return _beta3; };

		void calculateThresholdValues();
		double pre(uint nf) const;
		double post(uint nf) const;

		double evaluate(double Q0, double Qf, double alpha0) const;

		void update(uint nf);

		inline bool resumTabulated() const { return (_nf == _nff); }
		inline bool resumThreshold() const { return !resumTabulated(); }

	private:
		void assertNf() const;
		void assertScheme() const;
		
		double postMatch(double alpha, uint nf);
		double preMatch(double alpha, uint nf);

		double calcBeta0(uint nf) const;
		double calcBeta1(uint nf) const;
		double calcBeta2(uint nf) const;
		double calcBeta3(uint nf) const;
	};

} // class AlphaS

#endif // __ALPHAS_HPP
