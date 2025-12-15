#ifndef __DISTRIBUTION_HPP
#define __DISTRIBUTION_HPP

#include "Candia-v2/Common.hpp"

#include <array>
#include <numbers>

namespace Candia2
{

	class Distribution
	{
	protected:
		double _Q0; //!< chosen initial energy to evaluate alpha_s at
		double _alpha0; //!< value of alpha_s at chosen initial energy
		uint _nfi; //!< initial number of massless flavors
		std::array<double,8> _masses; //!< chosen quark masses
	public:

		Distribution() = delete; //!< must provide arguments
		Distribution(
			const double Q0, const double alpha0, const uint nfi,
			std::array<double,8> const& masses) 
			: _Q0(Q0), _alpha0(alpha0), _nfi(nfi), _masses(masses) 
		{ }
		virtual ~Distribution() = default;

		inline double Q0() const { return _Q0; }
		inline double alpha0() const { return _alpha0; }
		inline uint nfi() const { return _nfi; }
		inline std::array<double,8> const& masses() const { return _masses; }
		inline double masses(const uint idx) const { return _masses[idx]; }

		virtual double xuv(const double x) const = 0;
		virtual double xdv(const double x) const = 0;
		virtual double xg (const double x) const = 0;
		virtual double xdb(const double x) const = 0;
		virtual double xub(const double x) const = 0;
		virtual double xs (const double x) const = 0;
	};


	class LesHouchesDistribution final : public Distribution
	{
	private:
		static constexpr std::array<double,8> _leshouche_masses =
			// x    u    d            s                    c            b     t     x                    
			{ 0.0, 0.0, 0.0, std::numbers::sqrt2, std::numbers::sqrt2, 4.5, 175.0, 0.0 };
			// 'x' is a placeholder in the above array, not Bjorken-x
	public:
		LesHouchesDistribution()
			: Distribution(std::numbers::sqrt2, 0.35, 3, _leshouche_masses)
		{ }
		
		double xuv(const double x) const override;
		double xdv(const double x) const override;
		double xg (const double x) const override;
		double xdb(const double x) const override;
		double xub(const double x) const override;
		double xs (const double x) const override;
	}; // class Distribution
	
} // namespace Candia2


#endif // __DISTRIBUTION_HPP
