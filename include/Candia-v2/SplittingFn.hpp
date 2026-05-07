/**
 *  @file SplittingFn.hpp
 *  @brief Contains the @a SplittingFunction class, a derivation of @a Expression, to handle the splitting functions.
 */

#pragma once

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"
#include "Candia-v2/SpecialFuncs.hpp"

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
		static double _log_muf2_mur2;    //!< log of mu_f/mu_r
		P3ApproxType _imod{P3ApproxType::ImodAvg}; //!< approximation type for n3lo splitting functions
		static double flagNLP; //!< random flag enclosed in exact p3 splitting functions
	public:
		SplittingFunction() = default; //!< default constructor
		virtual ~SplittingFunction() = default; //!< default deconstructor

		SplittingFunction(P3ApproxType imod) : _imod{imod}
		{
		}

		/** updates the global value of nf and \f$\beta_0\f$ */
		inline static void update(uint nf, double beta0, double log_muf2_mur2)
		{
			_nf = nf;
			_beta0 = beta0;
		    _log_muf2_mur2 = log_muf2_mur2;
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
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};

	/**
	 *  @ingroup lo_splitfuncs
	 *  @brief Implements \f$P_{qq}^{(0)}\f$
	 */
	class P0qq final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};

	/**
	 *  @ingroup lo_splitfuncs
	 *  @brief Implements \f$P_{qg}^{(0)}\f$
	 */
	class P0qg final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
	};

	/**
	 *  @ingroup lo_splitfuncs
	 *  @brief Implements \f$P_{gq}^{(0)}\f$
	 */
	class P0gq final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
	};

	/**
	 *  @ingroup lo_splitfuncs
	 *  @brief Implements \f$P_{gg}^{(0)}\f$
	 */
	class P0gg final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};


	/**
	 *  @ingroup nlo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{+,(1)}\f$
	 */
	class P1nsp final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};

	/**
	 *  @ingroup nlo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{-,(1)}\f$
	 */
	class P1nsm final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};

	/**
	 *  @ingroup nlo_splitfuncs
	 *  @brief Implements \f$P_{qq}^{(1)}\f$
	 */
	class P1qq final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};

	/**
	 *  @ingroup nlo_splitfuncs
	 *  @brief Implements \f$P_{qg}^{(1)}\f$
	 */
	class P1qg final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
	};

	/**
	 *  @ingroup nlo_splitfuncs
	 *  @brief Implements \f$P_{gq}^{(1)}\f$
	 */
	class P1gq final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
	};

	/**
	 *  @ingroup nlo_splitfuncs
	 *  @brief Implements \f$P_{gg}^{(1)}\f$
	 */
	class P1gg final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};




	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{+,(2)}\f$
	 */
	class P2nsp final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{-,(2)}\f$
	 */
	class P2nsm final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{V,(2)}\f$
	 */
	class P2nsv final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{ps}}^{(2)}\f$ (pure-singlet)
	 */
	class P2ps final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{qq}^{(2)}\f$
	 */
	class P2qq final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{qg}^{(2)}\f$
	 */
	class P2qg final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{gq}^{(2)}\f$
	 */
	class P2gq final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
	};

	/**
	 *  @ingroup nnlo_splitfuncs
	 *  @brief Implements \f$P_{gg}^{(2)}\f$
	 */
	class P2gg final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};



	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{+,(3)}\f$
	 */
	class P3nsp final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{-,(3)}\f$
	 */
	class P3nsm final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{NS}}^{V,(3)}\f$
	 */
	class P3nsv final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{\mathrm{ps}}^{(3)}\f$ (pure-singlet)
	 */
	class P3ps final: public SplittingFunction
	{
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{qq}^{(3)}\f$
	 */
	class P3qq final: public SplittingFunction
	{
	private:
		double x1l4cff, x1l3cff;
		double y1l4cff, y1l3cff;
		double bfkl1;
		double x0l6cff, x0l5cff, x0l4cff;
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
		void preCalc() override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{qg}^{(3)}\f$
	 */
	class P3qg final: public SplittingFunction
	{
	private:
		double x1l5cff, x1l4cff;
		double y1l5cff, y1l4cff;
		double bfkl1;
		double x0l6cff, x0l5cff, x0l4cff;
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		void preCalc() override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{gq}^{(3)}\f$
	 */
	class P3gq final: public SplittingFunction
	{
	private:
		double x1l5cff, x1l4cff;
		double y1l5cff, y1l4cff;
		double bfkl0, bfkl1;
		double x0l6cff, x0l5cff, x0l4cff;
	public:
		using SplittingFunction::SplittingFunction;

		double calcRegular(double x) const override;
		void preCalc() override;
	};

	/**
	 *  @ingroup n3lo_splitfuncs
	 *  @brief Implements \f$P_{gg}^{(3)}\f$
	 */
	class P3gg final: public SplittingFunction
	{
	private:
		double a4gluon;
		double ccoeff, dcoeff;
		double x1l4cff, x1l3cff;
		double bfkl0, bfkl1;
		double x0l6cff, x0l5cff, x0l4cff;
	public:
		using SplittingFunction::SplittingFunction;
		
		double calcRegular(double x) const override;
		double calcPlus() const override;
		double calcDelta() const override;
		void preCalc() override;
	};


	/** @defgroup nnlo_splitfuncs_fortran Fortran Versions of P2 Splitting Functions */
	namespace mvv_p2
	{
		/**
		 *  @ingroup nnlo_splitfuncs_fortran
		 *  @brief Implements \f$P_{\mathrm{NS}}^{+,(2)}\f$
		 */
		class P2nsp final : public SplittingFunction
		{
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				return p2nspa(&x, &nf)/8.0;
			}
			inline double calcPlus() const override
			{
				double x = 0.0;
				int nf = static_cast<int>(_nf);
				return p2nsb(&x, &nf)/8.0;
			}
			inline double calcDelta() const override
			{
				double x = 1.0;
				int nf = static_cast<int>(_nf);
				return p2nspc(&x, &nf)/8.0;
			}
		};

		/**
		 *  @ingroup nnlo_splitfuncs_fortran
		 *  @brief Implements \f$P_{\mathrm{NS}}^{-,(2)}\f$
		 */
		class P2nsm final : public SplittingFunction
		{
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				return p2nsma(&x, &nf)/8.0;
			}
			inline double calcPlus() const override
			{
				double x = 0.0;
				int nf = static_cast<int>(_nf);
				return p2nsb(&x, &nf)/8.0;
			}
			inline double calcDelta() const override
			{
				double x = 1.0;
				int nf = static_cast<int>(_nf);
				return p2nsmc(&x, &nf)/8.0;
			}
		};

		/**
		 *  @ingroup nnlo_splitfuncs_fortran
		 *  @brief Implements \f$P_{\mathrm{NS}}^{V,(2)}\f$
		 */
		class P2nsv final : public SplittingFunction
		{
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				return (p2nsma(&x, &nf) + p2nssa(&x, &nf))/8.0;
			}
			inline double calcPlus() const override
			{
				double x = 0.0;
				int nf = static_cast<int>(_nf);
				return p2nsb(&x, &nf)/8.0;
			}
			inline double calcDelta() const override
			{
				double x = 1.0;
				int nf = static_cast<int>(_nf);
				return p2nsmc(&x, &nf)/8.0;
			}
		};

		/**
		 *  @ingroup nnlo_splitfuncs_fortran
		 *  @brief Implements \f$P_{qq}^{(2)}\f$
		 */
		class P2qq final : public SplittingFunction
		{
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				return (p2nspa(&x, &nf) + p2psa(&x, &nf))/8.0;
			}
			inline double calcPlus() const override
			{
				double x = 0.0;
				int nf = static_cast<int>(_nf);
				return p2nsb(&x, &nf)/8.0;
			}
			inline double calcDelta() const override
			{
				double x = 1.0;
				int nf = static_cast<int>(_nf);
				return p2nspc(&x, &nf)/8.0;
			}
		};

		/**
		 *  @ingroup nnlo_splitfuncs_fortran
		 *  @brief Implements \f$P_{qg}^{(2)}\f$
		 */
		class P2qg final : public SplittingFunction
		{
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				return p2qga(&x, &nf)/8.0;
			}
		};

		/**
		 *  @ingroup nnlo_splitfuncs_fortran
		 *  @brief Implements \f$P_{gq}^{(2)}\f$
		 */
		class P2gq final : public SplittingFunction
		{
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				return p2gqa(&x, &nf)/8.0;
			}
		};

		/**
		 *  @ingroup nnlo_splitfuncs_fortran
		 *  @brief Implements \f$P_{gg}^{(2)}\f$
		 */
		class P2gg final : public SplittingFunction
		{
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				return p2gga(&x, &nf)/8.0;
			}
			inline double calcPlus() const override
			{
				double x = 0.0;
				int nf = static_cast<int>(_nf);
				return p2ggb(&x, &nf)/8.0;
			}
			inline double calcDelta() const override
			{
				double x = 1.0;
				int nf = static_cast<int>(_nf);
				return p2ggc(&x, &nf)/8.0;
			}
		};
	} // namespace mvv_p2

	/** @defgroup n3lo_splitfuncs_exact Fortran Versions of P3 Splitting Functions */
	namespace mvv_p3
	{
		/**
		 *  @ingroup n3lo_splitfuncs_fortran
		 *  @brief Implements \f$P_{\mathrm{NS}}^{+,(3)}\f$
		 */
		class P3nsp final : public SplittingFunction
		{
		public:
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3nspa(&x, &nf, &imod)/16.0;
			}
			inline double calcPlus() const override
			{
				double x = 0.0;
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3nspb(&x, &nf, &imod)/16.0;
			}
			inline double calcDelta() const override
			{
				double x = 1.0;
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3nspc(&x, &nf, &imod)/16.0;
			}
		};

		/**
		 *  @ingroup n3lo_splitfuncs_fortran
		 *  @brief Implements \f$P_{\mathrm{NS}}^{-,(3)}\f$
		 */
		class P3nsm final: public SplittingFunction
		{
		public:
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3nsma(&x, &nf, &imod)/16.0;
			}
			inline double calcPlus() const override
			{
				double x = 0.0;
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3nsmb(&x, &nf, &imod)/16.0;
			}
			inline double calcDelta() const override
			{
				double x = 1.0;
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3nsmc(&x, &nf, &imod)/16.0;
			}
		};

		/**
		 *  @ingroup n3lo_splitfuncs_fortran
		 *  @brief Implements \f$P_{\mathrm{NS}}^{V,(3)}\f$
		 */
		class P3nsv final: public SplittingFunction
		{
		public:
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return (p3nsma(&x, &nf, &imod) + p3nssa(&x, &nf, &imod))/16.0;
			}
			inline double calcPlus() const override
			{
				double x = 0.0;
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3nsmb(&x, &nf, &imod)/16.0;
			}
			inline double calcDelta() const override
			{
				double x = 1.0;
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3nsmc(&x, &nf, &imod)/16.0;
			}
		};

		/**
		 *  @ingroup n3lo_splitfuncs_fortran
		 *  @brief Implements \f$P_{\mathrm{ps}}^{(3)}\f$ (pure-singlet)
		 */
		class P3ps final: public SplittingFunction
		{
		public:
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3psa(&x, &nf, &imod)/16.0;
			}
		};

		/**
		 *  @ingroup n3lo_splitfuncs_fortran
		 *  @brief Implements \f$P_{qq}^{(3)}\f$
		 */
		class P3qq final: public SplittingFunction
		{
		public:
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return (p3nspa(&x, &nf, &imod) + p3psa(&x, &nf, &imod))/16.0;
			}
			inline double calcPlus() const override
			{
				double x = 0.0;
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3nspb(&x, &nf, &imod)/16.0;
			}
			inline double calcDelta() const override
			{
				double x = 1.0;
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3nspc(&x, &nf, &imod)/16.0;
			}
		};

		/**
		 *  @ingroup n3lo_splitfuncs_fortran
		 *  @brief Implements \f$P_{qg}^{(3)}\f$
		 */
		class P3qg final: public SplittingFunction
		{
		public:
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3qga(&x, &nf, &imod)/16.0;
			}
		};

		/**
		 *  @ingroup n3lo_splitfuncs_fortran
		 *  @brief Implements \f$P_{gq}^{(3)}\f$
		 */
		class P3gq final: public SplittingFunction
		{
		public:
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3gqa_2512(&x, &nf, &imod)/16.0;
			}
		};

		/**
		 *  @ingroup n3lo_splitfuncs_fortran
		 *  @brief Implements \f$P_{gg}^{(3)}\f$
		 */
		class P3gg final: public SplittingFunction
		{
		public:
			using SplittingFunction::SplittingFunction;
		
			inline double calcRegular(double x) const override
			{
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3gga_2410(&x, &nf, &imod)/16.0;
			}
			inline double calcPlus() const override
			{
				double x = 0.0;
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3ggb_2410(&x, &nf, &imod)/16.0;
			}
			inline double calcDelta() const override
			{
				double x = 1.0;
				int nf = static_cast<int>(_nf);
				int imod = static_cast<int>(_imod);
				return p3ggc_2410(&x, &nf, &imod)/16.0;
			}
		};
	} // namespace mvv_p3

	/** @defgroup n3lo_splitfuncs_fortran Fortran Versions of P3 Splitting Functions */
	namespace p3_exact
	{
		/**
		 *  @ingroup n3lo_splitfuncs_exact
		 *  @brief Implements \f$P_{\mathrm{NS}}^{+,(3)}\f$
		 */
		class P3nsp final: public SplittingFunction
		{
		public:
			using SplittingFunction::SplittingFunction;
		
			double calcRegular(double x) const override;
			double calcPlus() const override;
			double calcDelta() const override;
		};

		/**
		 *  @ingroup n3lo_splitfuncs_exact
		 *  @brief Implements \f$P_{\mathrm{NS}}^{-,(3)}\f$
		 */
		class P3nsm final: public SplittingFunction
		{
		public:
			using SplittingFunction::SplittingFunction;
		
			double calcRegular(double x) const override;
			double calcPlus() const override;
			double calcDelta() const override;
		};

		/**
		 *  @ingroup n3lo_splitfuncs_exact
		 *  @brief Implements \f$P_{\mathrm{NS}}^{V,(3)}\f$
		 */
		class P3nsv final: public SplittingFunction
		{
		public:
			using SplittingFunction::SplittingFunction;
		
			double calcRegular(double x) const override;
			double calcPlus() const override;
			double calcDelta() const override;
		};
	}
	
} // namespace Candia2

