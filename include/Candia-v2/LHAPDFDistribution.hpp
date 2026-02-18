#ifndef __LHAPDFDISTRIBUTION_HPP
#define __LHAPDFDISTRIBUTION_HPP

#include "Candia-v2/Distribution.hpp"

#include <memory>

#include <LHAPDF/LHAPDF.h>

namespace Candia2
{
	using lhapdf_pdf_ptr_type = std::unique_ptr<LHAPDF::PDF>; //!< alias for a unique_ptr of an LHAPDF::PDF
	/** Helper function to return an @a lhapdf_pdf_ptr_type from standard LHAPDF::PDF arguments. */
	static inline lhapdf_pdf_ptr_type make_lhapdf_pdf(std::string const& setname, int num=0)
	{
		return lhapdf_pdf_ptr_type(LHAPDF::mkPDF(setname, num));
	}

	/**
	 *  @brief Class that implements initial conditions that are taken from an LHAPDF PDF.
	 */
	class LHAPDFDistribution final : public Distribution
	{
	public:
		using typename Distribution::value_type;
		using typename Distribution::accessor_type;
	private:
		lhapdf_pdf_ptr_type _pdf; //!< underlying PDF type
		std::vector<int> _pids; //!< list of the pdf's PIDs
	public:
		/**
		 *  @brief constructs a standard LHAPDFDistribution
		 *  @param lhapdf_pdf underlying pdf to use
		 *  @param Q0 the initial evolution energy
		 *  @param Qf the final evolution energy
		 */
		LHAPDFDistribution(lhapdf_pdf_ptr_type lhapdf_pdf, value_type Q0, value_type Qf);

		/** Getter for underlying pdf */
		inline LHAPDF::PDF const& pdf() const { return *_pdf; }

		value_type xg (value_type x) const override { return _pdf->xfxQ(21, x, _Q0); }
		value_type xu (value_type x) const override { return _pdf->xfxQ(1,  x, _Q0); }
		value_type xd (value_type x) const override { return _pdf->xfxQ(2,  x, _Q0); }
		value_type xs (value_type x) const override { return _pdf->xfxQ(3,  x, _Q0); }
		value_type xc (value_type x) const override { return _pdf->xfxQ(4,  x, _Q0); }
		value_type xdb(value_type x) const override { return _pdf->xfxQ(-1, x, _Q0); }
		value_type xub(value_type x) const override { return _pdf->xfxQ(-2, x, _Q0); }
		value_type xsb(value_type x) const override { return _pdf->xfxQ(-3, x, _Q0); }
		value_type xcb(value_type x) const override { return _pdf->xfxQ(-4, x, _Q0); }

		void fillSingletCoeffs(
			accessor_type const& accessor,
			std::vector<value_type> const& grid_points) const override;
		void fillNonSingletCoeffs(
			accessor_type const& accessor,
			std::vector<value_type> const& grid_points) const override;
	};
}

#endif // __LHAPDFDISTRIBUTION_HPP
