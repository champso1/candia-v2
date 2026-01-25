#ifndef __DISTRIBUTION_HPP
#define __DISTRIBUTION_HPP

#include "Candia-v2/Common.hpp"

#include <array>
#include <memory>
#include <numbers>

#include "LHAPDF/LHAPDF.h"

namespace Candia2
{

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

		virtual double xuv(double x) const = 0;
		virtual double xdv(double x) const = 0;
		virtual double xg (double x) const = 0;
		virtual double xdb(double x) const = 0;
		virtual double xub(double x) const = 0;
		virtual double xs (double x) const = 0;

		inline virtual double xu(double x) const { return xuv(x) + xub(x); }
		inline virtual double xd(double x) const { return xdv(x) + xdb(x); }
		inline virtual double xqplus(double x) const
		{ 
			return xuv(x) + 2.0*xub(x)
				+ xdv(x) + 2.0*xdb(x)
				+ 2.0*xs(x);
		}
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
		{}

		double xuv(double x) const override;
		double xdv(double x) const override;
		double xg (double x) const override;
		double xdb(double x) const override;
		double xub(double x) const override;
		double xs (double x) const override;
	}; 

	// TODO: this might be unnecessary
	using lhapdf_pdf_ptr_type = std::unique_ptr<LHAPDF::PDF>;

	class LHAPDFDistribution final : public Distribution
	{
	private:
		lhapdf_pdf_ptr_type _pdf; //!< actual lhapdf pdf object
		std::vector<int> _pids;

	public:
		LHAPDFDistribution(lhapdf_pdf_ptr_type lhapdf_pdf, double Q0);

		static inline lhapdf_pdf_ptr_type make_lhapdf_pdf(std::string const& setname, int num=0)
		{
			return lhapdf_pdf_ptr_type(LHAPDF::mkPDF(setname, num));
		}

		inline LHAPDF::PDF const& pdf() const { return *_pdf; }

		double xuv(double x) const override;
		double xdv(double x) const override;
		double xg (double x) const override;
		double xdb(double x) const override;
		double xub(double x) const override;
		double xs (double x) const override;
	};
	
} // namespace Candia2


#endif // __DISTRIBUTION_HPP
