/**
 *  @file OperatorMatrixElements.hpp
 *  @brief Contains the @a OpMatElem class, a derivation of @a Expression, to handle the operator matrix elements.
 */

#ifndef __OPERATOR_MATRIX_ELEMENTS_HPP
#define __OPERATOR_MATRIX_ELEMENTS_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"

#include <ome/ome.h>

namespace Candia2
{
	/**
	 *  @brief Class to handle the implementation of the operator matrix elements.
	 */
	class OpMatElem : public Expression
	{
	protected:
		static double _lm; //!< log(m_h^2/mu_r^2) = -log_mur2_muf2  ** NOTE THE MINUS **
		static uint _nf;   //!< number of active/massless flavors

		OpMatElem() = default; //!< default constructor
	public:
		virtual ~OpMatElem() = default; //!< default destructor

		/** Updates the internally stored value of @a lm and @a nf */
		inline static void update(double lm, uint nf)
		{
			log(LOG_INFO, "OME", "Setting L_M = {}, nf = {}", lm, nf);
			_lm = lm;
			_nf = nf;
		}
	};

	/**
	 *  @defgroup nnloopmatelems NNLO Operator Matrix Elements
	 *  @defgroup n3loopmatelems N3LO Operator Matrix Elements
	 */

	/**
	 *  @brief Implements \f$A_{qq,h}^{\mathrm{NS},(2)}\f$
	 *  @ingroup nnloopmatelems
	 */
	class A2ns final : public OpMatElem
	{
	public:
	    double calcRegular(double x) const override;
		double calcPlus(double x) const override;
		double calcDelta(double x) const override;
	};

	/**
	 *  @brief Implements \f$A_{gq,h}^{(2)}\f$
	 *  @ingroup nnloopmatelems
	 */
	class A2gq final : public OpMatElem
	{
	public:
		double calcRegular(double x) const override;
	};

	/**
	 *  @brief Implements \f$A_{gg,h}^{(2)}\f$
	 *  @ingroup nnloopmatelems
	 */
	class A2gg final : public OpMatElem
	{
	public:
		double calcRegular(double x) const override;
		double calcPlus(double x) const override;
		double calcDelta(double x) const override;
	};

	/**
	 *  @brief Implements \f$A_{hq,h}^{(2)}\f$
	 *  @ingroup nnloopmatelems
	 */
	class A2hq final : public OpMatElem
	{
	public:
		double calcRegular(double x) const override;
	};

	/**
	 *  @brief Implements \f$A_{hg,h}^{(2)}\f$
	 *  @ingroup nnloopmatelems
	 */
	class A2hg final : public OpMatElem
	{
	public:
		double calcRegular(double x) const override;
	};

	

	/**
	 *  @brief Implements the N3LO operator matrix elements using an underlying libome interface.
	 *  @ingroup n3loopmatelems
	 */
	class OpMatElemN3LO final : public OpMatElem
	{
	public:
		/** Alias for underlying libome type */
		using ome_type = ome::rpd_distribution<ome::ome_as_view<double>, ome::ome_as_plus_view<double>, ome::ome_as_const_view<double>>;

	private:
		ome_type _ome; //!< underlying libome interface

	public:
		/** constructs an @a OpMatElemN3LO with the specified libome type */
		OpMatElemN3LO(ome_type const& ome) : _ome{ome} {}
		~OpMatElemN3LO() = default; //!< default deconstructor

		inline double calcRegular(double x) const override
		{
			if (!_ome.has_regular())
				return 0;

			auto reg = _ome.get_regular().value();
			return reg[3](_lm, _nf, x);
		}

		inline double calcPlus(double x) const override
		{
			if (!_ome.has_plus())
				return 0;

			auto plus = _ome.get_plus().value();
			return plus[3](_lm, _nf, x);
		}

		inline double calcDelta(double x) const override
		{
			UNUSED(x);

			if (!_ome.has_delta())
				return 0;

			auto delta = _ome.get_delta().value();
			return delta[3](_lm, _nf);
		}
	};

	extern "C"
	{
		double aqqpss3_(double* x, double* NF, double* lm);
	}

	class AQqPSs3 final : public OpMatElem
	{
	public:
		inline double calcRegular(double x) const override
		{
			double nf = static_cast<double>(_nf);
			return aqqpss3_(&x, &nf, &_lm);
		}
	};
};


#endif // __OPERATOR_MATRIX_ELEMENTS_HPP
