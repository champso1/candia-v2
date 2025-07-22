#ifndef __OPERATOR_MATRIX_ELEMENTS_HPP
#define __OPERATOR_MATRIX_ELEMENTS_HPP

#include "Candia-v2/Expression.hpp"

namespace Candia2
{

	/** Class to represent the transition operator matrix elements,
	 *  also known as the matching conditions
	 */
	class OpMatElem : public Expression
	{
	public:
		/** @name Constructors/destructors
		 */
		///@{
		OpMatElem() = default;
		~OpMatElem() = default;
		///@}
	};



	class A2ns : public OpMatElem
	{
	public:
		A2ns() = default;
		~A2ns()  = default;

		double Regular(const double x) const;
		double Plus(const double x) const;
		double Delta(const double x) const;
	};

	class A2gq : public OpMatElem
	{
	public:
		A2gq() = default;
		~A2gq()  = default;

		double Regular(const double x) const;
	};

	class A2gg : public OpMatElem
	{
	public:
		A2gg() = default;
		~A2gg()  = default;

		double Regular(const double x) const;
		double Plus(const double x) const;
		double Delta(const double x) const;
	};

	class A2hq : public OpMatElem
	{
	public:
		A2hq() = default;
		~A2hq()  = default;

		double Regular(const double x) const;
	};

	class A2hg : public OpMatElem
	{
	public:
		A2hg() = default;
		~A2hg()  = default;

		double Regular(const double x) const;
	};
	
};


#endif // __OPERATOR_MATRIX_ELEMENTS_HPP
