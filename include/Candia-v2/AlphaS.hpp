/**
 *  @file AlphaS.hpp
 *  @brief Contains the @a AlphaS class which handles the evolution of the QCD running coupling.
 */

#ifndef __ALPHAS_HPP
#define __ALPHAS_HPP

#include <array>

#include "Candia-v2/Common.hpp"

namespace Candia2
{

	class LHAPDFDistribution;

	/**
	 *  @brief Class to handle the evolution of the QCD running coupling
	 */
	class AlphaS
	{
	private:
		uint _order{}; //!< Perturbative order.
		uint _nf{};    //!< number of active flavors (set via setFFNS())
		uint _nfi{}, _nff{}; //!< initial and final number of flavors
		double _Q0{}, _Qf{}; //!< initial and final evolution energies
		double _alpha0{}; //!< Initial value of \f$\alpha_s\f$ at Q0

	    double _beta0{}; //!< \f$\beta_0\f$
		double _beta1{}; //!< \f$\beta_1\f$
		double _beta2{}; //!< \f$\beta_2\f$
		double _beta3{}; //!< \f$\beta_3\f$

		double _L{}; //!< log (mu_R/mu_F)

		std::array<double, 8> _masses{}; //!< Values of quark masses.

		std::array<double, 8> _pre{};  //!< array of pre-threshold values of \f$\alpha_s\f$
		std::array<double, 8> _post{}; //!< array of post-threshold values of \f$\alpha_s\f$

		/** Enum for determing if we are in VFNS or FFNS */
	    enum Scheme : int
		{
			UNSET = -1,
			FIXED = 0,
			VARIABLE = 1
		};
		int _scheme; //!< current scheme
	public:
		/**
		 *  @brief Initializes initial conditions for \f$\alpha_s\f$.
		 *  @param order The perturbative order.
		 *  @param Q0 The initial evolution energy.
		 *  @param Qf The final evolution energy.
		 *  @param alpha0 The value of \f$\alpha_s\f$ at the initial energy.
		 *  @param log_mur2_muf2 The ratio of the renormalization scale \f$\mu_r^2\f$ to the factorization scale \f$\mu_f^2\f$.
		 */
	    AlphaS(uint order, double Q0, double Qf, double alpha0, double log_mur2_muf2)
			: _order{order}, _Q0{Q0}, _Qf{Qf}, _alpha0{alpha0}, _L{log_mur2_muf2}, _scheme{UNSET}
		{}

		/**
		 *  @brief Sets \f$\alpha_s\f$ to evolve with a fixed number of flavors \f$n_f\f$
		 *  @param nf The number of flavors to fix.
		 */
		void setFFNS(uint nf);

		/**
		 *  @brief Sets \f$\alpha_s\f$ to evolve with a variable number of flavors, with thresholds determined by @a masses.
		 *  @param masses The array of on-shell masses at which \f$\alpha_s\f$ will perform matching
		 *  @param nfi The starting value of \f$n_f\f$, determined by the selected initial distribution.
		 */
		void setVFNS(std::array<double, 8> const& masses, uint nfi);

		double masses(uint nf) const; //!< getter for the mass corresponding to @a nf
		inline uint nfi() const { return _nfi; } //!< getter for the starting value of \f$n_f\f$
		inline uint nff() const { return _nff; } //!< getter for final value of \f$n_f\f$, determined by final evolution energy
	    
		double betaFn(double alpha) const; //!< returns the value of the beta-function given \f$\alpha\f$ as @a alpha

		/** @brief retrieves the stored value of \f$\beta_0\f$ */
		inline double beta0() const { return _beta0; };
		/** @brief retrieves the stored value of \f$\beta_1\f$ */
		inline double beta1() const { return _beta1; };
		/** @brief retrieves the stored value of \f$\beta_2\f$ */
		inline double beta2() const { return _beta2; };
		/** @brief retrieves the stored value of \f$\beta_3\f$ */
		inline double beta3() const { return _beta3; };

		void calculateThresholdValues(); //!< given the mass array, calculates the value of \f$\alpha-s\f$ pre and post-threshold
		double pre(uint nf) const; //!< returns the value of \f$\alpha_s\f$ before the threshold corresponding to \f$n_f\f$
		double post(uint nf) const; //!< returns the value of \f$\alpha_s\f$ after the threshold corresponding to \f$n_f\f$

		/** @brief evaluates \f$\alpha_s\f$ using the solution to the Cauchy problem.
		 *
		 *  Has an exact solution at LO, otherwise performs a fourth-order runge-kutta method
		 *
		 *  @param Q0 the initial value at which \f$\alpha_s\f$ is evaluated at
		 *  @param Qf the final value to evaluate \f$\alpha_s\f$ at
		 *  @param alpha0 the value of \f$\alpha_s\f$ at the initial energy @a Q0
		 */
		double evaluate(double Q0, double Qf, double alpha0) const;

		/**
		 *  @brief Updates the values of the beta-coefficients at @a nf
		 *  @param nf requested value of \f$n_f\f$
		 */
		void update(uint nf);

		/** @brief Returns true if we are resumming to a tabluated energy, false otherwise */
		inline bool resumTabulated() const { return (_nf == _nff); }
		/** @brief Returns true if we are resumming to a threshold energy, false otherwise */
		inline bool resumThreshold() const { return !resumTabulated(); }

	private:
		void assertNf() const; //!< asserts whether the current value of @a nf is valid
		void assertScheme() const;  //!< asserts whether the user's choice(s) given the current scheme is(are) valid
		// 
		double postMatch(double alpha, uint nf); //!< calculates \f$\alpha_s\f$ post threshold given the value pre-threshold
		double preMatch(double alpha, uint nf); //!< calculates \f$\alpha_s\f$ pre threshold given the value post-threshold

		double calcBeta0(uint nf) const; //!< calculates \f$\beta_0\f$ for @a nf active flavors
		double calcBeta1(uint nf) const; //!< calculates \f$\beta_1\f$ for @a nf active flavors
		double calcBeta2(uint nf) const; //!< calculates \f$\beta_2\f$ for @a nf active flavors
		double calcBeta3(uint nf) const; //!< calculates \f$\beta_3\f$ for @a nf active flavors
	};

} // class AlphaS

#endif // __ALPHAS_HPP
