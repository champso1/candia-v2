/** @file
 *
 *  All splitting function information is stored here.
 *  Every splitting function contains a "Regular" piece,
 *  a singlar piece related to the plus distribution (called "Plus")
 *  and a piece related to the delta function (called "Delta")
 */

#ifndef __SPLITTINGFN_HPP
#define __SPLITTINGFN_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"

namespace Candia2
{	
	/** Class to represent a generic splitting function/kernel.
	 */
	class SplittingFunction : public Expression
	{
	protected:
		static uint _nf; //!< number of active/currently massless flavors
		static double _beta0; //!< beta0 coefficient (for P0gg)
		static double _kr; //!< log of mu_f/mu_r

		SplittingFunction() = default;

		// static HPLog _hplog; //!< hplog interface object
	public:
		virtual ~SplittingFunction() = default;

		/** @name Setters/getters
		 */
		///@{
		inline static uint Nf() { return _nf; }
		inline static double Beta0() { return _beta0; }
		inline static double kr() { return _kr; }
		inline static void kr(double kr) { _kr = kr; }
		///@}
		
		
		/** @brief Updates the stored value of nf and associated beta0
		 */
		inline static void Update(const uint nf, const double beta0)
		{ _nf = nf; _beta0 = beta0; }
		

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
	class P0ns final : public SplittingFunction
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
	class P0qq final: public SplittingFunction
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
	class P0qg final: public SplittingFunction
	{
	public:
		P0qg() { }
		~P0qg() = default;

		double Regular(const double x) const override;
	};

	/** Leading order, g->q singlet kernel.
	 */
	class P0gq final: public SplittingFunction
	{
	public:
		P0gq() { }
		~P0gq() = default;

		double Regular(const double x) const override;
	};

	/** Leading order, g->g singlet kernel.
	 */
	class P0gg final: public SplittingFunction
	{
	public:
		P0gg() { }
		~P0gg() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};

	///@}

	/** @name Next-to leading order (NLO) splitting functions
	 */
	///@{

	/** Next-to leading order, non-singlet kernel (plus component)
	 */
	class P1nsp final: public SplittingFunction
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
	class P1nsm final: public SplittingFunction
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
	class P1qq final: public SplittingFunction
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
	class P1qg final: public SplittingFunction
	{
	public:
		P1qg() { }
		~P1qg() = default;

		double Regular(const double x) const override;
	};

	/** Next-to leading order, g->q singlet kernel.
	 */
	class P1gq final: public SplittingFunction
	{
	public:
		P1gq() { }
		~P1gq() = default;

		double Regular(const double x) const override;
	};

	/** Next-to leading order, g->g singlet kernel.
	 */
	class P1gg final: public SplittingFunction
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
	class P2nsp final: public SplittingFunction
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
	class P2nsm final: public SplittingFunction
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
	class P2nsv final: public SplittingFunction
	{
	public:
		P2nsv() { }
		~P2nsv() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};


	class P2ps final: public SplittingFunction
	{
	public:
		P2ps() = default;
		~P2ps() = default;

		double Regular(const double x) const override;
		// double Plus(const double x) const override;
		// double Delta(const double x) const override;
	};




	/** Next-to-next-to leading order, q->q singlet kernel.
	 */
	class P2qq final: public SplittingFunction
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
	class P2qg final: public SplittingFunction
	{
	public:
		P2qg() { }
		~P2qg() = default;

		double Regular(const double x) const override;
	};

	/** Next-to-next-to leading order, g->q singlet kernel.
	 */
	class P2gq final: public SplittingFunction
	{
	public:
		P2gq() { }
		~P2gq() = default;

		double Regular(const double x) const override;
	};

	/** Next-to-next-to leading order, g->g singlet kernel.
	 */
	class P2gg final: public SplittingFunction
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

	/** N3LO non-singlet kernel (plus component)
	 */
	class P3nsp final: public SplittingFunction
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

	/** N3LO non-singlet kernel (minus component)
	 */
	class P3nsm final: public SplittingFunction
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

	
	/** N3LO non-singlet kernel (valence component)
	 */
	class P3nsv final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use
		
	public:
		P3nsv(const uint imod=3) : _imod(imod) { }
		~P3nsv() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};


	/** N3LO pure-singlet kernel
	 */
	class P3ps final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use

	public:
		P3ps(const uint imod=3) : _imod(imod) { }
		~P3ps() = default;

		double Regular(const double x) const override;
		// double Plus(const double x) const override;
		// double Delta(const double x) const override;
	};

	/** N3LO qq diagonal singlet term (P3nsp + P3ps)
	 */ 
	class P3qq final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use

	public:
		P3qq(const uint imod=3) : _imod(imod) { }
		~P3qq() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};



	/** N3LO qg off-diagonal singlet term
	 */
	class P3qg final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use

	public:
		P3qg(const uint imod=3) : _imod(imod) { }
		~P3qg() = default;

		double Regular(const double x) const override;
		//double Plus(const double x) const override;
		//double Delta(const double x) const override;
	};

	/** N3LO gq off-diagonal singlet term
	 */
	class P3gq final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use

	public:
		P3gq(const uint imod=3) : _imod(imod) { }
		~P3gq() = default;

		double Regular(const double x) const override;
		// double Plus(const double x) const override;
		// double Delta(const double x) const override;
	};

	/** N3LO gg diagonal singlet term
	 */
	class P3gg final: public SplittingFunction
	{
	private:
		const uint _imod; //!< flag for which approximation to use

	public:
		P3gg(const uint imod=3) : _imod(imod) { }
		~P3gg() = default;

		double Regular(const double x) const override;
		double Plus(const double x) const override;
		double Delta(const double x) const override;
	};	
	///@}







	namespace splitfunc_fortran
	{

		/** @name Next-to-next-to leading order (NNLO) splitting functions.
		 */
		///@{

		/** Next-to-next-to leading order, non-singlet kernel (plus component)
		 */
		class P2nsp final: public SplittingFunction
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
		class P2nsm final: public SplittingFunction
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
		class P2nsv final: public SplittingFunction
		{
		public:
			P2nsv() { }
			~P2nsv() = default;

			double Regular(const double x) const override;
			double Plus(const double x) const override;
			double Delta(const double x) const override;
		};


		class P2ps final: public SplittingFunction
		{
		public:
			P2ps() = default;
			~P2ps() = default;

			double Regular(const double x) const override;
			// double Plus(const double x) const override;
			// double Delta(const double x) const override;
		};




		/** Next-to-next-to leading order, q->q singlet kernel.
		 */
		class P2qq final: public SplittingFunction
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
		class P2qg final: public SplittingFunction
		{
		public:
			P2qg() { }
			~P2qg() = default;

			double Regular(const double x) const override;
		};

		/** Next-to-next-to leading order, g->q singlet kernel.
		 */
		class P2gq final: public SplittingFunction
		{
		public:
			P2gq() { }
			~P2gq() = default;

			double Regular(const double x) const override;
		};

		/** Next-to-next-to leading order, g->g singlet kernel.
		 */
		class P2gg final: public SplittingFunction
		{
		public:
			P2gg() { }
			~P2gg() = default;

			double Regular(const double x) const override;
			double Plus(const double x) const override;
			double Delta(const double x) const override;
		};
	}

} // namespace Candia2


#endif // __SPLITTINGFN_HPP
