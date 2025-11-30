#ifndef __OPERATOR_MATRIX_ELEMENTS_HPP
#define __OPERATOR_MATRIX_ELEMENTS_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"
#include "ome/ome_type_aliases.h"
#include "ome/rpd_distribution.h"

#include <iostream>

#include <ome/ome.h>

namespace Candia2
{

	/** Class to represent the transition operator matrix elements,
	 *  also known as the matching conditions
	 */
	class OpMatElem : public Expression
	{
	protected:
		static double _lm; //!< log(m_h^2/mu_r^2) = -log_mur2_muf2  ** NOTE THE MINUS **
		static uint _nf;   //!< number of flavors

		OpMatElem() = default;
	public:
		virtual ~OpMatElem() = default;

		/** @brief updates the values of @a lm and @a nf
		 */
		inline static void update(double lm, uint nf)
		{
			std::cerr << "[OME] Setting L_M = " << lm << ", nf = " << nf << '\n';
			_lm = lm;
			_nf = nf;
		}
	};



	class A2ns final : public OpMatElem
	{
	public:
		A2ns() = default;
		~A2ns()  = default;

		double Regular(double x) const override;
		double Plus(double x) const override;
		double Delta(double x) const override;
	};

	class A2gq final : public OpMatElem
	{
	public:
		A2gq() = default;
		~A2gq()  = default;

		double Regular(double x) const override;
	};

	class A2gg final : public OpMatElem
	{
	public:
		A2gg() = default;
		~A2gg()  = default;

		double Regular(double x) const override;
		double Plus(double x) const override;
		double Delta(double x) const override;
	};

	class A2hq final : public OpMatElem
	{
	public:
		A2hq() = default;
		~A2hq()  = default;

		double Regular(double x) const override;
	};

	class A2hg final : public OpMatElem
	{
	public:
		A2hg() = default;
		~A2hg()  = default;

		double Regular(double x) const override;
	};


	class OpMatElemN3LO final : public OpMatElem
	{
		using ome_type = ome::rpd_distribution<ome::ome_as_view<double>, ome::ome_as_plus_view<double>, ome::ome_as_const_view<double>>;
		ome_type _ome;

	public:
		OpMatElemN3LO(ome_type const& ome) : _ome{ome} {}
		~OpMatElemN3LO() = default;

		inline double Regular(double x) const override
		{
			if (!_ome.has_regular())
				return 0;

			auto reg = _ome.get_regular().value();
			return reg[3](_lm, _nf, x);
		}

		inline double Plus(double x) const override
		{
			if (!_ome.has_plus())
				return 0;

			auto plus = _ome.get_plus().value();
			return plus[3](_lm, _nf, x);
		}

		inline double Delta(double x) const override
		{
			UNUSED(x);

			if (!_ome.has_delta())
				return 0;

			auto delta = _ome.get_delta().value();
			return delta[3](_lm, _nf);
		}
	};
};


#endif // __OPERATOR_MATRIX_ELEMENTS_HPP
