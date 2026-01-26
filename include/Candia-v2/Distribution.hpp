#ifndef __DISTRIBUTION_HPP
#define __DISTRIBUTION_HPP

#include "Candia-v2/Common.hpp"

#include <array>
#include <memory>
#include <numbers>

#include "LHAPDF/LHAPDF.h"

namespace Candia2
{
    using coefficient_accessor_type = std::function<double&(uint,uint)>;

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
		Distribution() = default;	
		Distribution(double Q0, double a0, uint nfi, masses_type const& masses)	
			: _Q0{Q0}, _alpha0{a0}, _nfi{nfi}, _masses{masses}
		{}
		virtual ~Distribution() = default;

		inline double Q0() const { return _Q0; }
		inline double alpha0() const { return _alpha0; }
		inline uint nfi() const { return _nfi; }
		inline masses_type const& masses() const { return _masses; }
		inline double masses(uint idx) const { return _masses[idx]; }

		// initial conditions
		// the defaults in terms of virtuality are in accordance with the standard benchmark initial conditions
		// in which the gluon, u/d/s (and anti) have initial conditions, while the others dont
		// the charm is possible since it has a smallish mass, but by default is just 0
		// the top is of course not even included
		virtual double xg(double x)  const = 0;
		virtual double xu(double x)  const = 0;
		virtual double xd(double x)  const = 0;
		virtual double xs(double x)  const = 0;
		virtual double xc(double x)  const { return 0.0; }
		virtual double xub(double x) const = 0;
		virtual double xdb(double x) const = 0;
		virtual double xsb(double x) const = 0;
		virtual double xcb(double x) const { return 0.0; }

		virtual void fillSingletCoeffs(
			coefficient_accessor_type const& accessor,
			std::vector<double> const& grid_points) const = 0;
		virtual void fillNonSingletCoeffs(
			coefficient_accessor_type const& accessor,
			std::vector<double> const& grid_points) const = 0;

		inline virtual double xqplus(double x) const
		{ 
			return xu(x) + xub(x)
				+ xd(x) + xdb(x)
				+ xs(x) + xsb(x)
				+ xc(x) + xcb(x);
		}
	};


	class LesHouchesDistribution final : public Distribution
	{
	private:
		static constexpr std::array<double,8> _leshouche_masses = {	
		  // x    u    d            s                    c            b     t     x                    
			0.0, 0.0, 0.0, std::numbers::sqrt2, std::numbers::sqrt2, 4.5, 175.0, 0.0 };
		// 'x' is a placeholder in the above array, not Bjorken-x
	public:
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

	using lhapdf_pdf_ptr_type = std::unique_ptr<LHAPDF::PDF>;
	static inline lhapdf_pdf_ptr_type make_lhapdf_pdf(std::string const& setname, int num=0)
	{
		return lhapdf_pdf_ptr_type(LHAPDF::mkPDF(setname, num));
	}

	class LHAPDFDistribution final : public Distribution
	{
	private:
		lhapdf_pdf_ptr_type _pdf;
		std::vector<int> _pids;

	public:
		LHAPDFDistribution(lhapdf_pdf_ptr_type lhapdf_pdf, double Q0, double Qf);

		inline LHAPDF::PDF const& pdf() const { return *_pdf; }

		double xg (double x) const override { return _pdf->xfxQ(21, x, _Q0); }
		double xu (double x) const override { return _pdf->xfxQ(1,  x, _Q0); }
		double xd (double x) const override { return _pdf->xfxQ(2,  x, _Q0); }
		double xs (double x) const override { return _pdf->xfxQ(3,  x, _Q0); }
		double xc (double x) const override { return _pdf->xfxQ(4,  x, _Q0); }
		double xdb(double x) const override { return _pdf->xfxQ(-1, x, _Q0); }
		double xub(double x) const override { return _pdf->xfxQ(-2, x, _Q0); }
		double xsb(double x) const override { return _pdf->xfxQ(-3, x, _Q0); }
		double xcb(double x) const override { return _pdf->xfxQ(-4, x, _Q0); }

		void fillSingletCoeffs(
			coefficient_accessor_type const& accessor,
			std::vector<double> const& grid_points) const override;
		void fillNonSingletCoeffs(
			coefficient_accessor_type const& accessor,
			std::vector<double> const& grid_points) const override;
		
	};
	
} // namespace Candia2


#endif // __DISTRIBUTION_HPP
