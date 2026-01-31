/**
 *  @file SplittingFn.hpp
 *  @brief Contains the @a SplittingFunction class, a derivation of @a Expression, to handle the splitting functions.
 */

#ifndef __SPLITTINGFN_HPP
#define __SPLITTINGFN_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"

namespace Candia2
{
	/**
	 *  @brief Class to handle the implementation of the splitting functions.
	 */
	class SplittingFunction : public Expression
	{
	protected:
		static uint _nf;      //!< number of active/currently massless flavors
		static double _beta0; //!< \f$\beta_0\f$ value for P0gg calculation
		static double _kr;    //!< log of mu_f/mu_r
	public:
		SplittingFunction() = default; //!< default constructor
		virtual ~SplittingFunction() = default; //!< default deconstructor

		/** updates the global value of nf and \f$\beta_0\f$ */
		inline static void update(uint nf, double beta0)
		{
			_nf = nf; _beta0 = beta0;
		}
	};


	/** @defgroup lo_splitfuncs Leading-order (LO) splitting functions */
	/** @defgroup nlo_splitfuncs Next-to-leading-order (NLO) splitting functions */
	/** @defgroup nnlo_splitfuncs Next-to-next-to-leading-order splitting functions */
	/** @defgroup n3lo_splitfuncs Next-to-next-to-next-to-leading-order (N3LO) splitting functions */

	/**
	 *  @ingroup lo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{(0)}\f$
	 */
	class P0ns final : public SplittingFunction
	{
	public:
	    using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup lo_splitfuncs
	 *  @brief Implements \f$P_{qq}^{(0)}\f$
	 */
	class P0qq final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup lo_splitfuncs
	 *  @brief Implements \f$P_{qg}^{(0)}\f$
	 */
	class P0qg final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
	};

	/**
	 *  @ingroup lo_splitfuncs
	 *  @brief Implements \f$P_{gq}^{(0)}\f$
	 */
	class P0gq final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
	};

	/**
	 *  @ingroup lo_splitfuncs
	 *  @brief Implements \f$P_{gg}^{(0)}\f$
	 */
	class P0gg final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};


	/**
	 *  @ingroup nlo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{+,(1)}\f$
	 */
	class P1nsp final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup nlo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{-,(1)}\f$
	 */
	class P1nsm final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup nlo_splitfuncs
	 *  @brief Implements \f$P_{qq}^{(1)}\f$
	 */
	class P1qq final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup nlo_splitfuncs
	 *  @brief Implements \f$P_{qg}^{(1)}\f$
	 */
	class P1qg final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
	};

	/**
	 *  @ingroup nlo_splitfuncs
	 *  @brief Implements \f$P_{gq}^{(1)}\f$
	 */
	class P1gq final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
	};

	/**
	 *  @ingroup nlo_splitfuncs
	 *  @brief Implements \f$P_{gg}^{(1)}\f$
	 */
	class P1gg final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};




	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{+,(2)}\f$
	 */
	class P2nsp final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{-,(2)}\f$
	 */
	class P2nsm final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{V,(2)}\f$
	 */
	class P2nsv final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{ps}}^{(2)}\f$ (pure-singlet)
	 */
	class P2ps final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{qq}^{(2)}\f$
	 */
	class P2qq final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{qg}^{(2)}\f$
	 */
	class P2qg final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{gq}^{(2)}\f$
	 */
	class P2gq final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{gg}^{(2)}\f$
	 */
	class P2gg final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};



	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{+,(3)}\f$
	 */
	class P3nsp final: public SplittingFunction
	{
	private:
		const uint _imod{3}; //!< flag for which approximation to use
		
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{-,(3)}\f$
	 */
	class P3nsm final: public SplittingFunction
	{
	private:
		const uint _imod{3}; //!< flag for which approximation to use
		
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{V,(3)}\f$
	 */
	class P3nsv final: public SplittingFunction
	{
	private:
		const uint _imod{3}; //!< flag for which approximation to use
		
	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{ps}}^{(3)}\f$ (pure-singlet)
	 */
	class P3ps final: public SplittingFunction
	{
	private:
		const uint _imod{3}; //!< flag for which approximation to use

	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{qq}^{(3)}\f$
	 */
	class P3qq final: public SplittingFunction
	{
	private:
		const uint _imod{3}; //!< flag for which approximation to use

	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{qg}^{(3)}\f$
	 */
	class P3qg final: public SplittingFunction
	{
	private:
		const uint _imod{3}; //!< flag for which approximation to use

	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{gq}^{(3)}\f$
	 */
	class P3gq final: public SplittingFunction
	{
	private:
		const uint _imod{3}; //!< flag for which approximation to use

	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{gg}^{(3)}\f$
	 */
	class P3gg final: public SplittingFunction
	{
	private:
		const uint _imod{3}; //!< flag for which approximation to use

	public:
		using SplittingFunction::SplittingFunction;
		
		double _reg_func(double x) const override;
		double _plus_func(double x) const override;
		double _delta_func(double x) const override;
	};
} // namespace Candia2


#endif // __SPLITTINGFN_HPP
