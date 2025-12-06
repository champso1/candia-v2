/** @file
 *
 *  The main file containing the main evolution code. This interface is what the user
 *  will interface with 99% of the time.
 */

#ifndef __CANDIA_HPP
#define __CANDIA_HPP


#include "Candia-v2/Common.hpp"
#include "Candia-v2/Grid.hpp"
#include "Candia-v2/AlphaS.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/SplittingFn.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"
#include "Candia-v2/FuncArrGrid.hpp"

#include <concepts>
#include <cstdlib>
#include <functional>
#include <memory>
#include <complex>
#include <utility>
#include <vector>

namespace Candia2
{	
	/** @brief Main entry point for numerically solving the DGLAP equations
	 */
	class DGLAPSolver
	{
	private:
		uint _order{}; //!< perturbative order
		Grid _grid; //!< main @a Grid object
		double _Qf{}; //!< final energy value to evolve to
		AlphaS _alpha_s; //!< main @a AlphaS object
		double _mur2_muf2{}; //!< mu_r^2/mu_f^2
		double _log_mur2_muf2{}; //!< log(_mur2_muf2)

		uint _qct{}; //!< counter for tabulated Q's (to be removed)
		
		std::unique_ptr<Distribution> _dist; //!< main (initial) distribution

		/** @name Active flavor counts
		 */
		///@{
		uint _nf{}; //!< current active number of massless flavors
		uint _nfi{}; //!< minimum number based on initial evolution and provided quark masses
		uint _nff{}; //!< final based on final evolution and provided quark masses
		///@}

		/** @name Alpha_s threshold values
		 *  values of alpha post-match of the previous interval
		 *  and pre-match of the next interval,
		 *  i.e. \f$\alpha_s\f$ in the beginning/end of the current interval
		 */
		///@{
		double _alpha0{};
		double _alpha1{};
		///@}

		/** @name other evolution variables
		 */
		///@{
		uint _iterations{}; //!< number of singlet/non-singlet iterations
		uint _trunc_idx{}; //!< number of additional singlet truncated iterations
		///@}

		
		/** @name non-singlet coefficients
		 */
		///@{
	    MultiDimVector<double, 3>::type _A{}; //!< LO Coefficients 
		MultiDimVector<double, 4>::type _B{}; //!< NLO Coefficients
		MultiDimVector<double, 5>::type _C{}; //!< NNLO Coefficients
		MultiDimVector<double, 6>::type _D{}; //!< N3LO Coefficients

		MultiDimArrayGrid_t<2> _A2{};
		MultiDimArrayGrid_t<3> _B2{};
		MultiDimArrayGrid_t<4> _C2{};
		MultiDimArrayGrid_t<5> _D2{};
		///@}

		/** @name singlet coefficients
		 */
		///@{
		MultiDimVector<double, 4>::type _S{}; //!< singlet coefficients
	    MultiDimVector<ArrayGrid, 3>::type _S2{};
		///@}

		MultiDimVector<double, 2>::type _F{}; //!< final distribution
		std::vector<ArrayGrid> _F2{}; //!< final distribution ALT


		/** @name N3LO additional parameters
		 */
		///@{
		std::array<double,8> _r1{}; //!< one real solution to N3LO quadratic
		std::array<std::complex<double>,8> _r2{};
		std::array<std::complex<double>,8> _r3{};
		std::array<double,8> _b{};  //!< -2*Re[r2]
		std::array<double,8> _c{};  //!< |r2|^2
		///@}

		/** @brief list of all splitting functions
		 */
		std::map<std::string, FunctionGrid> _expression_grids{};

		template <typename TExpr, typename... TExprArgs>
		requires std::derived_from<TExpr, Expression>
		void createExpressionGrid(std::string const& name, Grid const& grid, TExprArgs&&... args)
		{
			_expression_grids.emplace(
				std::make_pair(
					name, FunctionGrid(grid, std::unique_ptr<Expression>(new TExpr(std::forward<TExprArgs>(args)...)))
				)
			);
		}
		
		inline
	    FunctionGrid&
		getSplitFunc(std::string const& name)
		{
			auto splitfunc = _expression_grids.find(name);
			if (splitfunc == _expression_grids.end())
			{
				std::cerr << "[DGLAP2] Splitting function '" << name << "' does not exist.\n";
				exit(EXIT_FAILURE);
			}
			return _expression_grids.at(name);
		}

		inline 
		FunctionGrid&
		getOME(std::string const& name)
		{
			return _expression_grids.at(name);
		}

		void loadAllExpressions();
		
		/** @defgroup SplitFuncs Splitting function pointers
		 */
		///@{
		///@}

		/** @name LO splitting functions
		 *  @ingroup SplitFuncs
		 */
		///@{
		std::shared_ptr<SplittingFunction> _P0ns{};
		std::shared_ptr<SplittingFunction> _P0qq{};
		std::shared_ptr<SplittingFunction> _P0qg{};
		std::shared_ptr<SplittingFunction> _P0gq{};
		std::shared_ptr<SplittingFunction> _P0gg{};
		///@}

		/** @name NLO splitting functions
		 *  @ingroup SplitFuncs
		 */
		///@{
		std::shared_ptr<SplittingFunction> _P1nsp{};
		std::shared_ptr<SplittingFunction> _P1nsm{};
		std::shared_ptr<SplittingFunction> _P1qq{};
		std::shared_ptr<SplittingFunction> _P1qg{};
		std::shared_ptr<SplittingFunction> _P1gq{};
		std::shared_ptr<SplittingFunction> _P1gg{};
		///@}

		/** @name NNLO splitting functions
		 *  @ingroup SplitFuncs
		 */
		///@{
		std::shared_ptr<SplittingFunction> _P2nsp{};
		std::shared_ptr<SplittingFunction> _P2nsm{};
		std::shared_ptr<SplittingFunction> _P2nsv{};
		std::shared_ptr<SplittingFunction> _P2qq{};
		std::shared_ptr<SplittingFunction> _P2qg{};
		std::shared_ptr<SplittingFunction> _P2gq{};
		std::shared_ptr<SplittingFunction> _P2gg{};
		///@}

		/** @name NNNLO splitting functions
		 *  @ingroup SplitFuncs
		 */
		///@{
		std::shared_ptr<SplittingFunction> _P3nsp{};
		std::shared_ptr<SplittingFunction> _P3nsm{};
		std::shared_ptr<SplittingFunction> _P3nsv{};
		std::shared_ptr<SplittingFunction> _P3ps{};
		std::shared_ptr<SplittingFunction> _P3qq{};
		std::shared_ptr<SplittingFunction> _P3qg{};
		std::shared_ptr<SplittingFunction> _P3gq{};
		std::shared_ptr<SplittingFunction> _P3gg{};
		///@}


		/** @name Operator Matrix Elements
		 */
		///@{
		std::shared_ptr<OpMatElem> _A2ns{};
		std::shared_ptr<OpMatElem> _A2gq{};
		std::shared_ptr<OpMatElem> _A2gg{};
		std::shared_ptr<OpMatElem> _A2hq{};
		std::shared_ptr<OpMatElem> _A2hg{};

		std::shared_ptr<OpMatElem> _A3nsp{};
		std::shared_ptr<OpMatElem> _A3nsm{};
		std::shared_ptr<OpMatElem> _A3gq{};
		std::shared_ptr<OpMatElem> _A3gg{};
		std::shared_ptr<OpMatElem> _A3hq{};
		std::shared_ptr<OpMatElem> _A3hg{};
		std::shared_ptr<OpMatElem> _A3psqq{};
		std::shared_ptr<OpMatElem> _A3sqg{};
		///@}


		/** @defgroup RecRels Recursion relations
		 */
		///@{
		///@}

		/** @name Singlet recursion relations
		 *  @ingroup RecRels
		 */
		///@{
		double RecRelS_1(
			std::vector<double> const& S,
			uint k,
			std::shared_ptr<SplittingFunction> P);
		double RecRelS_1(
			ArrayGrid & S,
			uint k,
			FunctionGrid & P);

		double RecRelS_2(
			std::vector<double> const& S1,
			std::vector<double> const& S2,
			uint k,
			std::shared_ptr<SplittingFunction> P0,
			std::shared_ptr<SplittingFunction> P1);
		double RecRelS_2(
			ArrayGrid & S1,
			ArrayGrid & S2,
			uint k,
			FunctionGrid & P0,
			FunctionGrid & P1);

		double RecRelS_3(
			std::vector<double> const& S1,
			std::vector<double> const& S2,
			std::vector<double> const& S3,
			uint k,
			std::shared_ptr<SplittingFunction> P0,
			std::shared_ptr<SplittingFunction> P1,
			std::shared_ptr<SplittingFunction> P2);
		double RecRelS_3(
			ArrayGrid & S1,
			ArrayGrid & S2,
			ArrayGrid & S3,
			uint k,
			FunctionGrid & P0,
			FunctionGrid & P1,
			FunctionGrid & P2);

		double RecRelS_4(
			std::vector<double> const& S1,
			std::vector<double> const& S2,
			std::vector<double> const& S3,
			std::vector<double> const& S4,
			uint k,
			std::shared_ptr<SplittingFunction> P0,
			std::shared_ptr<SplittingFunction> P1,
			std::shared_ptr<SplittingFunction> P2,
			std::shared_ptr<SplittingFunction> P3);
		double RecRelS_4(
			ArrayGrid & S1,
			ArrayGrid & S2,
			ArrayGrid & S3,
			ArrayGrid & S4,
			uint k,
			FunctionGrid & P0,
			FunctionGrid & P1,
			FunctionGrid & P2,
			FunctionGrid & P3);
		///@}

		/** @name LO recursion relations
		 *  @ingroup RecRels
		 */
		///@{
		double RecRelLO(
			std::vector<double> const& A,
			uint k,
			std::shared_ptr<SplittingFunction> P0);

		double RecRelLO(
			ArrayGrid & A,
			uint k,
			FunctionGrid & P0);
		///@}

		/** @name NLO recursion relations
		 *  @ingroup RecRels
		 */
		///@{
		double RecRelNLO_1(
			std::vector<double> const& b,
			uint k,
			std::shared_ptr<SplittingFunction> P0);
		double RecRelNLO_2(
			std::vector<double> const& B,
			uint k,
			std::shared_ptr<SplittingFunction> P0,
			std::shared_ptr<SplittingFunction> P1);

		double RecRelNLO_1(
			ArrayGrid & b,
			uint k,
			FunctionGrid & P0);
		double RecRelNLO_2(
			ArrayGrid & B,
			uint k,
			FunctionGrid & P0,
			FunctionGrid & P1);
		///@}

		/** @name NNLO recursion relations
		 *  @ingroup RecRels
		 */
		///@{
		double RecRelNNLO_1(
			std::vector<double> const& C,
			uint k,
			std::shared_ptr<SplittingFunction> P0);
		double RecRelNNLO_2(
			std::vector<double> const& C,
			uint k,
			std::shared_ptr<SplittingFunction> P0,
			std::shared_ptr<SplittingFunction> P1,
			std::shared_ptr<SplittingFunction> P2);
		double RecRelNNLO_3(
			std::vector<double> const& C,
			uint k,
			std::shared_ptr<SplittingFunction> P0,
			std::shared_ptr<SplittingFunction> P1);

		double RecRelNNLO_1(
			ArrayGrid & C,
			uint k,
			FunctionGrid& P0);
		double RecRelNNLO_2(
			ArrayGrid & C,
			uint k,
			FunctionGrid& P0,
			FunctionGrid& P1,
			FunctionGrid& P2);
		double RecRelNNLO_3(
			ArrayGrid & C,
			uint k,
			FunctionGrid& P0,
			FunctionGrid& P1);
		///@}

		/** @name NNNLO recursion relations
		 *  @ingroup RecRels
		 */
		///@{
		double RecRelN3LO_1(
			std::vector<double> const& D,
			uint k,
			std::shared_ptr<SplittingFunction> P0);
		double RecRelN3LO_2(
			std::vector<double> const& D,
			uint k,
			std::shared_ptr<SplittingFunction> P0,
			std::shared_ptr<SplittingFunction> P1,
			std::shared_ptr<SplittingFunction> P2,
			std::shared_ptr<SplittingFunction> P3);
		double RecRelN3LO_3(
			std::vector<double> const& D,
			uint k,
			std::shared_ptr<SplittingFunction> P0,
			std::shared_ptr<SplittingFunction> P1,
			std::shared_ptr<SplittingFunction> P2,
			std::shared_ptr<SplittingFunction> P3);
		double RecRelN3LO_4(
			std::vector<double> const& D,
			uint k,
			std::shared_ptr<SplittingFunction> P0,
			std::shared_ptr<SplittingFunction> P1,
			std::shared_ptr<SplittingFunction> P2,
			std::shared_ptr<SplittingFunction> P3);

		double RecRelN3LO_1(
			ArrayGrid& D,
			uint k,
			FunctionGrid& P0);
		double RecRelN3LO_2(
			ArrayGrid& D,
			uint k,
			FunctionGrid& P0,
			FunctionGrid& P1,
			FunctionGrid& P2,
			FunctionGrid& P3);
		double RecRelN3LO_3(
			ArrayGrid& D,
			uint k,
			FunctionGrid& P0,
			FunctionGrid& P1,
			FunctionGrid& P2,
			FunctionGrid& P3);
		double RecRelN3LO_4(
			ArrayGrid& D,
			uint k,
			FunctionGrid& P0,
			FunctionGrid& P1,
			FunctionGrid& P2,
			FunctionGrid& P3);
		///@}

		bool _multi_thread;
	public:
		/** @name Constructors/destructors
		 */
		///@{
		/** @brief Main constructor to initialize evolution
			@param order: perturbative order
		 *  @param grid: initialized grid object
		 *  @param Qf: final energy to evaluate to
		 *  @param iterations: number of iterations
		 *  @param trunc_idx: number of additional truncated iterations (singlet) per main iteration
		 *  @param initial_dist: initial distribution for masses, pdfs, and alpha0
		 *  @param kr: ratio of mu_r to mu_f
		 */
	    DGLAPSolver(const uint order, Grid & grid, const double Qf,
			const uint iterations, const uint trunc_idx,
			std::unique_ptr<Distribution> initial_dist,
			double kr = 1.0, bool multi_thread = false
		);
		~DGLAPSolver();
		///@}
		
		/** @brief Returns a const-ref to the \f$\alpha_s\f$ object.
		 */
		inline AlphaS const& GetAlphaS() const { return _alpha_s; }

		/** @brief Returns a const-ref to the @a Grid object.
		 */
		inline Grid const& GetGrid() const { return _grid; }

		/** @brief sets the number of normal and truncation iterations
		 */
		void SetEvolutionVariables(const uint iterations, const uint trunc_idx);

		
		/** @brief Main function to evolve coefficients.
		 *  @return The final distributions
		 */
		MultiDimVector<double, 2>::type Evolve();

		auto Evolve2() -> decltype(_F2);

		inline
	    auto const&
		getF2() const
		{
			return _F2;
		}

	private:
		/** @brief initializes coefficient arrays with distributions
		 */
		void SetInitialConditions();

		/** @brief For NNLO and beyond, evaluates the splitting functions at all grid points
		 *  and caches the results for speed
		 */
		void FillSplittingFunctionCaches();
		
		/** @brief Computes the actual distributions that are evolved.
		 */
		void SetupCoefficients();
		void SetupCoefficients2();

		void FixDistributions2(bool resum_tab, bool resum_threshold, std::vector<ArrayGrid>& temp_arr, std::vector<ArrayGrid>& temp_arr_singlet);

		/** @brief Computes (from the calculated coefficients)
		 */
		[[maybe_unused]]
		void SetupFinalDistributions();


		/** @brief evolves the singlet coefficients
		 */
		void EvolveSinglet();

		/** @brief evolves the non-singlet coefficients
		 */
		void EvolveNonSinglet();

		/** @brief evolves the non-singlet coefficients
		 *  @note does this with multi-threading, and makes
		 *  a new thread for each distribution
		 */
		void EvolveNonSingletThreaded();

		/** @brief evolves the singlet coefficients (different form of the code including n3lo)
		 */
		void EvolveSingletAlt();

		void EvolveSinglet2(std::reference_wrapper<std::vector<ArrayGrid>> arr, double L1);
		void EvolveNonSinglet2(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			double L1, double L2, double L3, double L4);
		void EvolveNonSinglet2Threaded(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			double L1, double L2, double L3, double L4);

		/** @brief After computing coefficients, do resummation to a tabulated energy value
		 *  @param Q: the tabulated energy value to evolve to
		 *  @return the final distributions
		 */
		MultiDimVector<double, 2>::type Resum(const double Q);

		/** @brief Finish resumming to the next threshold
		 */
		void ResumThreshold();

		/** Handles the matching conditions/transition operator matrix elements
		 */
		void HeavyFlavorTreatment();

		/** Handles the matching conditions/transition operator matrix elements ALTERNATE
		 */
		void HeavyFlavorTreatmentAlt();

		/** @brief helper functions for heavy flavor treatment
		 */
		///@{
		void HFT_NNLO1_ORIG(uint j, uint k);
		void HFT_NNLO2_ORIG(uint k);
		void HFT_NNLO3_ORIG(uint k);

		void HFT_N3LO1_ORIG(uint j, uint k);
		void HFT_N3LO2_ORIG(uint k);
		void HFT_N3LO3_ORIG(uint k);
		
		void HFT_NNLO1(uint j, uint k); //!< q_NS^+
		void HFT_NNLO2(uint j, uint k); //!< q_NS^-
		void HFT_NNLO3(uint k); //!< g
		void HFT_NNLO4(uint k); //!< (q_h + qb_h)
		void HFT_NNLO5(uint k); //!< q^+

		void HFT_N3LO1(uint j, uint k, double SP); //!< q_NS^+
		void HFT_N3LO2(uint j, uint k, double SP); //!< q_NS^-
		void HFT_N3LO3(uint k); //!< g
		void HFT_N3LO4(uint k); //!< (q_h + qb_h)
		void HFT_N3LO5(uint k); //!< q^+// 
		///@}

		void HeavyFlavorTreatment2();
		void HFT2_NNLO1(ArrayGrid& c, uint j, uint k);
		void HFT2_NNLO2(ArrayGrid& g, ArrayGrid& qp, uint k);
		void HFT2_NNLO3(ArrayGrid& g, ArrayGrid& qp, uint k);
		void HFT2_N3LO1(ArrayGrid& q, ArrayGrid& qb, uint j, uint k, double SP);
		void HFT2_N3LO2(ArrayGrid& q, ArrayGrid& qb, uint j, uint k, double SP);
		void HFT2_N3LO3(ArrayGrid& g, ArrayGrid& qp, uint k);
		void HFT2_N3LO4(ArrayGrid& g, ArrayGrid& qp, uint k);


	    void _mt_EvolveDistribution_NS_LO(uint j);
	    void _mt_EvolveDistribution_NS_NLO(uint j, std::shared_ptr<SplittingFunction> P1);
	    void _mt_EvolveDistribution_NS_NNLO(uint j, std::array<std::shared_ptr<SplittingFunction>, 2> const& P2);
	    void _mt_EvolveDistribution_NS_N3LO(uint j, std::array<std::shared_ptr<SplittingFunction>, 3> const& P3);

		void _mt2_EvolveDistribution_NS_LO  (
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, double L1);
	    void _mt2_EvolveDistribution_NS_NLO (
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::string const& P1, std::array<double, 2> const& L);
	    void _mt2_EvolveDistribution_NS_NNLO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::array<std::string, 2> const& P, std::array<double, 3> const& L);
	    void _mt2_EvolveDistribution_NS_N3LO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::array<std::string, 3> const& P, std::array<double, 4> const& L);

	public:
		enum LatexNumberFormat : uint
		{
			SCIENTIFIC,
			PERCENT
		};
		
		/** @brief Prints all required distributions into a nice latex table
		 */
		static void OutputLatexTable(MultiDimVector<double, 2>::type const& dists, std::string const& title, uint format);
	};

	

	
}


#endif // __CANDIA_HPP
