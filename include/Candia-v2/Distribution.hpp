/**
 *  @file Distribution.hpp
 *  @brief Contains the @a Distribution class and its derivations that implement various initial evolution conditions
 */

#ifndef __DISTRIBUTION_HPP
#define __DISTRIBUTION_HPP

#include "Candia-v2/Common.hpp"

#include <array>
#include <cmath>
#include <numbers>

namespace Candia2
{
	/**
	 *  A function object that takes a distribution and grid index
	 *  and returns a reference to a coefficient on the grid.
	 */
    using coefficient_accessor_type = std::function<double&(uint,uint)>;

	/**
	 *  @brief Base class for handling initial evolution conditions
	 */
	class Distribution
	{
	public:
		using masses_type = std::array<double, 8>;
	protected:
		double _Q0{}; //!< chosen initial energy to evaluate alpha_s at
		double _alpha0{}; //!< value of alpha_s at chosen initial energy
		uint _nfi{}; //!< initial number of massless flavors
		masses_type _masses{}; //!< chosen quark masses
	public:
		Distribution() = default; //!< default constructor for dists that manually setup values
		/**
		 *  @brief sets up some initial values for the distribution
		 *  @param Q0 the initial evolution energy
		 *  @param alpha0 the value of \f$\alpha_s\f$ at the initial evolution energy
		 *  @param nfi the number of initial massless flavors this distribution supports
		 *  @param masses array of on-shell quark masses
		 */
		Distribution(double Q0, double alpha0, uint nfi, masses_type const& masses)	
			: _Q0{Q0}, _alpha0{alpha0}, _nfi{nfi}, _masses{masses}
		{}
		virtual ~Distribution() = default;

		/** Getter for @a Q0 */
		inline double Q0() const { return _Q0; }
		/** Getter for @a alpha0 */
		inline double alpha0() const { return _alpha0; }
		/** Getter for @a nfi */
		inline uint nfi() const { return _nfi; }
		/** Getter for @a masses */
		inline masses_type const& masses() const { return _masses; }
		/** Getter for a mass in @a masses */
		inline double masses(uint idx) const { return _masses[idx]; }

		// initial conditions
		// the defaults in terms of virtuality are in accordance with the standard benchmark initial conditions
		// in which the gluon, u/d/s (and anti) have initial conditions, while the others dont
		// the charm is possible since it has a smallish mass, but by default is just 0
		// the bottom and top of course are not even included
		/** Retrieves the initial value of \f$xg(x)\f$ at @a x */
		virtual double xg(double x)  const = 0;
		/** Retrieves the initial value of \f$xu(x)\f$ at @a x */
		virtual double xu(double x)  const = 0;
		/** Retrieves the initial value of \f$xd(x)\f$ at @a x */
		virtual double xd(double x)  const = 0;
		/** Retrieves the initial value of \f$xs(x)\f$ at @a x */
		virtual double xs(double x)  const = 0;
		/** Retrieves the initial value of \f$xc(x)\f$ at @a x */
		virtual double xc(double x)  const { return 0.0; }
		/** Retrieves the initial value of \f$x\bar{u}(x)\f$ at @a x */
		virtual double xub(double x) const = 0;
		/** Retrieves the initial value of \f$x\bar{d}(x)\f$ at @a x */
		virtual double xdb(double x) const = 0;
		/** Retrieves the initial value of \f$x\bar{s}(x)\f$ at @a x */
		virtual double xsb(double x) const = 0;
		/** Retrieves the initial value of \f$x\bar{c}(x)\f$ at @a x */
		virtual double xcb(double x) const { return 0.0; }

		/**
		 *  @brief Helper for filling the set of singlet coefficients
		 *  @param accessor the accessor to retrieve a singlet coefficient
		 *  @param grid_points the grid points at which to fill the coefficients
		 */
		virtual void fillSingletCoeffs(
			coefficient_accessor_type const& accessor,
			std::vector<double> const& grid_points) const = 0;
		/**
		 *  @brief Helper for filling the set of non-singlet coefficients
		 *  @param accessor the accessor to retrieve a non-singlet coefficient
		 *  @param grid_points the grid points at which to fill the coefficients
		 */
		virtual void fillNonSingletCoeffs(
			coefficient_accessor_type const& accessor,
			std::vector<double> const& grid_points) const = 0;

		/** Returns the \f$xq^{(+)}\f$ distribution */
		inline virtual double xqplus(double x) const
		{ 
			return xu(x) + xub(x)
				+ xd(x) + xdb(x)
				+ xs(x) + xsb(x)
				+ xc(x) + xcb(x);
		}
	};

	/**
	 *  @brief Class that implements the standard benchmark initial conditions
	 */
	class LesHouchesDistribution final : public Distribution
	{
	private:
		/** array of mass points for this distribution */
		static constexpr std::array<double,8> _leshouche_masses = {	
		  // x    u    d            s                    c            b     t     x                    
			0.0, 0.0, 0.0, std::numbers::sqrt2, std::numbers::sqrt2, 4.5, 175.0, 0.0 };
		// 'x' is a placeholder in the above array, not Bjorken-x
	public:
		/**
		 *  @brief initializes the base @a Distribution class with Q0=\f$\sqrt{2}\f$, alpha0=0.35, nfi=3, and masses= @a _leshouche_masses
		 */
		LesHouchesDistribution()
			: Distribution(std::numbers::sqrt2, 0.35, 3, _leshouche_masses)
		{}

		inline double xuv(double x) const
		{
			return 5.1072*std::pow(x, 0.8)*std::pow(1.0-x, 3.0);
		}
		inline double xdv(double x) const
		{
			return 3.06432*std::pow(x, 0.8)*std::pow(1.0-x, 4.0);
		}

		inline double xg (double x) const override
		{
			return 1.7*std::pow(x, -0.1)*std::pow(1.0-x, 5.0);
		}
		inline double xu (double x) const override
		{
			return xuv(x) + xub(x);
		}
		inline double xd (double x) const override
		{
			return xdv(x) + xdb(x);
		}
		inline double xdb(double x) const override
		{
			return 0.1939875*std::pow(x, -0.1)*std::pow(1.0-x, 6.0);
		}
		inline double xub(double x) const override
		{
			return xdb(x)*(1.0-x);
		}
		inline double xs (double x) const override
		{
			return 0.2*(xub(x) + xdb(x));
		}
		inline double xsb(double x) const override { return xs(x); }

		void fillSingletCoeffs(
			coefficient_accessor_type const& accessor,
			std::vector<double> const& grid_points) const override;
	    void fillNonSingletCoeffs(
			coefficient_accessor_type const& accessor,
			std::vector<double> const& grid_points) const override;
	};
	
} // namespace Candia2


#endif // __DISTRIBUTION_HPP
