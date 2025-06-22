/** @file
 *
 *  Contains information related to the strong coupling constant \f$\alpha_s\f$
 *  and the coefficients of the \f$\beta\f$-function.
 */

#ifndef __ALPHAS_HPP
#define __ALPHAS_HPP

#include <vector>
#include <array>

#include "Candia-v2/Common.hpp"


namespace Candia2
{

	/** Class for handling the calculation and matching conditions
	 *  for the QCD strong coupling constant.
	 */
	class AlphaS
	{
	private:
		uint _order; //!< Perturbative order.
		double _Q0, _alpha0; //!< Initial value of \f$\alpha_s\f$ at Q0

		/** @name Stored values of beta coefficients/values
		 */
		///@{
	    double _beta0, _beta1, _beta2, _beta3;
		double _beta;
		///@}

		std::array<double, 8> _masses; //!< Values of quark masses.

		/** @name Values of \f$\alpha_s\f$ pre- and post-threshold.
		 */
		///@{
		std::array<double, 8> _pre, _post;
		///@}

	public:
		/** @brief Default constructor is deleted
		 */
		AlphaS() = delete;

		/** @brief AlphaS constructor
		 *  @param order: perturbative order
		 *  @param nf: number of massless flavors
		 *  @param Q0: initial factorization scale
		 *  @param alpha0: initial value of \f$\alpha_s\f$ at Q0
		 *  @param masses: list of quark masses
		 */
		AlphaS(const uint order, const uint nf,
			   const double Q0, const double alpha0,
			   std::array<double,8> const& masses)
			: _order(order), _Q0(Q0), _alpha0(alpha0), _masses(masses)
		{ Update(nf); }

		
		/** @brief getter for perturbative order
		 *  @return the perturbative order
		 */
		inline uint Order() const
		{ return _order; }
		
		/** @brief Getter for mass corresponding to flavor @a nf
		 */
		double Masses(const uint nf) const;
	    
		
		/** @brief Calculates the value of the beta function
		 *  @param alpha: current value of \f$\alpha_s\f$
		 *  @return The value of the beta function
		 */
		double BetaFn(const double alpha) const;


		/** @name The \f$\beta\f$-coefficients
		 *
		 *  These function simple return the already-calculated coeffs
		 *  These are updated by calling @a Update with a new
		 *  number of massless flavors
		 */
	    ///@{
		inline double Beta0() const { return _beta0; };
		inline double Beta1() const { return _beta1; };
		inline double Beta2() const { return _beta2; };
		inline double Beta3() const { return _beta3; };
		///@}
	    

		/**
		 */
		///@{
		/** @brief Calculates all threshold values of \f$\alpha_s\f$
		 *  @param Qf: the final energy/factorization scale
		 */
		void CalculateThresholdValues(const double Qf);
		
		
		/** @brief returns pre-match value for @nf
		 */
		double Pre(const uint nf) const;

		/** @brief returns post-match value for @nf
		 */
		double Post(const uint nf) const;
		///@}

		
		/** @brief Calculates the value of \f$\alpha_s\f$ at scale @a Qf given \f$\alpha_0\f$ at scale @a Q0
		 *
		 *  Utilizes an exact solution of the beta function at leading-order,
		 *  otherwise uses a fourth-order Runge-Kutta method.
		 *
		 *  @param Q0: the initial scale
		 *  @param Qf: the desired final scale
		 *  @param alpha0: the value of \f$\alpha_s\f$ at @a Q0
		 *  @return \f$\alpha_s\f$
		 */
		double Evaluate(const double Q0, const double Qf,
						const double alpha0) const;


		/** @brief Updates beta coefficients with new @a nf
		 */
		void Update(const uint nf);


	private:
		/** @name \f$\alpha_s\f$ threshold value calculations
		 */
		///@{
		/** @brief Calculates \f$\alpha_s\f$ in the next threshold
		 *  @param alpha: current \f$\alpha_s\f$
		 *  @return \f$\alpha_s\f$ in the next threshold
		 */
		double PostMatch(const double alpha);
		
		/** @brief Calculates \f$\alpha_s\f$ in the previous threshold
		 *  @param alpha: current \f$\alpha_s\f$
		 *  @return \f$\alpha_s\f$ in the previous threshold
		 */
		double PreMatch(const double alpha);

		
		/** @name The \f$\beta\f$ coefficients
		 */
	    ///@{
		double CalcBeta0(const uint nf) const;
		double CalcBeta1(const uint nf) const;
		double CalcBeta2(const uint nf) const;
		double CalcBeta3(const uint nf) const;
		////@}
	};

}

#endif // __ALPHAS_HPP
