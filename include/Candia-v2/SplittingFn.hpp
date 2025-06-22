/** @file
 *
 *  All splitting function information is stored here.
 *  Every splitting function contains a "Regular" piece,
 *  a singlar piece related to the plus distribution (called "Plus")
 *  and a piece related to the delta function (called "Delta")
 */

#ifndef __SPLITTINGFN_HPP
#define __SPLITTINGFN_HPP

#include <exception>
#include <memory>
#include <string>

#include "Candia-v2/Common.hpp"
#include "Candia-v2/HPLog.hpp"
#include "Candia-v2/Math.hpp"

namespace Candia2
{	
	/** Class to represent a generic splitting function/kernel.
	 */
	class SplittingFunction
	{
	protected:
		static uint _nf; //!< number of active/currently massless flavors

		// static HPLog _hplog; //!< hplog interface object
	public:
		/** @name Constructors/destructors
		 */
		///@{
		/** @brief Constructs a generic splitting function.
		 *  @param nf: Number of currently massless flavors (default 3)
		 */
		SplittingFunction() = default;
		~SplittingFunction() = default; //!< default destructor
		///@}

		/** @name Setters/getters
		 */
		///@{
		/** @brief Returns nf (a copy)
		 */
		inline uint Nf() const { return _nf; }

		/** @brief Returns a non-const reference to @a nf
		 */
		inline uint& Nf() { return _nf; }
		///@}

		/** @name Kernel components
		 */
		///@{
		/** @brief the regular part */
		virtual double Regular(const double x) const = 0;
		/** @brief the plus-distribution part */
		virtual double Plus(const double x) const = 0;
		/** @brief the delta-function part */
		virtual double Delta(const double x) const = 0;
		///@{


		

	protected:
		/** @name HPLOG interface functions
		 */
		///@{
		// double Li2(const double x) const;
		// double Li3(const double x) const;
		// double S12(const double x) const;
		///@}
	};


	/** @name Leading-order splitting functions.
	 */
	///@{

	/** Leading-order, non-singlet kernel.
	 */
	class P0ns : public SplittingFunction
	{
	public:
		P0ns() { }
		~P0ns() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};


	/** Leading order, q->q singlet kernel.
	 */
	class P0qq : public SplittingFunction
	{
	public:
		P0qq() { }
		~P0qq() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** Leading order, q->g singlet kernel.
	 */
	class P0qg : public SplittingFunction
	{
	public:
		P0qg() { }
		~P0qg() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** Leading order, g->q singlet kernel.
	 */
	class P0gq : public SplittingFunction
	{
	public:
		P0gq() { }
		~P0gq() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** Leading order, g->g singlet kernel.
	 */
	class P0gg : public SplittingFunction
	{
	public:
		P0gg(const double beta0) : _beta0(beta0) { }
		~P0gg() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;

	private:
		const double _beta0;
	};

	///@}

	/** @name Next-to leading order (NLO) splitting functions
	 */
	///@{

	/** Next-to leading order, non-singlet kernel (plus component)
	 */
	class P1nsp : public SplittingFunction
	{
	public:
		P1nsp() { }
		~P1nsp() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** Next-to leading order, non-singlet kernel (minus component)
	 */
	class P1nsm : public SplittingFunction
	{
	public:
		P1nsm() { }
		~P1nsm() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};


	/** Next-to leading order, q->q singlet kernel.
	 */
	class P1qq : public SplittingFunction
	{
	public:
		P1qq() { }
		~P1qq() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** Next-to leading order, q->g singlet kernel.
	 */
	class P1qg : public SplittingFunction
	{
	public:
		P1qg() { }
		~P1qg() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** Next-to leading order, g->q singlet kernel.
	 */
	class P1gq : public SplittingFunction
	{
	public:
		P1gq() { }
		~P1gq() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** Next-to leading order, g->g singlet kernel.
	 */
	class P1gg : public SplittingFunction
	{
	public:
		P1gg() { }
		~P1gg() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	///@}

	/** @name Next-to-next-to leading order (NNLO) splitting functions.
	 */
	///@{

	/** Next-to-next-to leading order, non-singlet kernel (plus component)
	 */
	class P2nsp : public SplittingFunction
	{
	public:
		P2nsp() { }
		~P2nsp() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** Next-to-next-to leading order, non-singlet kernel (minus component)
	 */
	class P2nsm : public SplittingFunction
	{
	public:
		P2nsm() { }
		~P2nsm() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** Next-to-next-to leading order, non-singlet kernel (valence component)
	 */
	class P2nsv : public SplittingFunction
	{
	public:
		P2nsv() { }
		~P2nsv() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};


	/** Next-to-next-to leading order, q->q singlet kernel.
	 */
	class P2qq : public SplittingFunction
	{
	public:
		P2qq() { }
		~P2qq() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** Next-to-next-to leading order, q->g singlet kernel.
	 */
	class P2qg : public SplittingFunction
	{
	public:
		P2qg() { }
		~P2qg() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** Next-to-next-to leading order, g->q singlet kernel.
	 */
	class P2gq : public SplittingFunction
	{
	public:
		P2gq() { }
		~P2gq() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** Next-to-next-to leading order, g->g singlet kernel.
	 */
	class P2gg : public SplittingFunction
	{
	public:
		P2gg() { }
		~P2gg() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	///@}


	/** @name Next-to-next-to-next-to (NNNLO or N3LO) splitting functions
	 */
	///@{

	/** NNNLO non-singlet kernel (plus component)
	 */
	class P3nsp : public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use
		
	public:
		P3nsp(const uint imod=3) : _imod(imod) { }
		~P3nsp() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** NNNLO non-singlet kernel (minus component)
	 */
	class P3nsm : public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use
		
	public:
		P3nsm(const uint imod=3) : _imod(imod)  { }
		~P3nsm() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	/** NNNLO non-singlet kernel (sea component)
	 */
	class P3nss : public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use
		
	public:
		P3nss(const uint imod=3) : _imod(imod) { }
		~P3nss() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};
	
	///@}

}; // namespace Candia2


#endif // __SPLITTINGFN_HPP
