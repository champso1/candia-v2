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

#include <chrono>
#include <memory>
#include <iostream>
#include <fstream>


namespace Candia2
{	
	/** @brief Main entry point for numerically solving the DGLAP equations
	 */
	class DGLAPSolver
	{
	private:
		uint _order; //!< perturbative order
		Grid _grid; //!< main @a Grid object
		double _Qf; //!< final energy value to evolve to
		AlphaS _alpha_s; //!< main @a AlphaS object

		uint _qct; //!< counter for tabulated Q's (to be removed)

		std::shared_ptr<Distribution> _dist; //!< main (initial) distribution

		/** @name Active flavor counts
		 */
		///@{
		uint _nf; //!< current active number of massless flavors
		uint _nfi; //!< minimum number based on initial evolution and provided quark masses
		uint _nff; //!< final based on final evolution and provided quark masses
		///@}

		/** @name Alpha_s threshold values
		 *  values of alpha post-match of the previous interval
		 *  and pre-match of the next interval,
		 *  i.e. \f$\alpha_s\f$ in the beginning/end of the current interval
		 */
		///@{
		double _alpha0;
		double _alpha1;
		///@}

		
		/** @name non-singlet coefficients
		 */
		///@{
	    MultiDimVector<double, 3>::type _A; //!< LO Coefficients
		MultiDimVector<double, 4>::type _B; //!< NLO Coefficients
		MultiDimVector<double, 5>::type _C; //!< NNLO Coefficients
		MultiDimVector<double, 6>::type _D; //!< N3LO Coefficients
		///@}

		/** @name singlet coefficients
		 */
		///@{
		MultiDimVector<double, 4>::type _S; //!< singlet coefficients
		///@}

		MultiDimVector<double, 2>::type _F; //!< final distribution


		/** @name N3LO additional parameters
		 */
		///@{
		std::array<double,8> _r1; //!< one real solution to N3LO quadratic
		std::array<double,8> _b; //!< -2*Re[r2]
		std::array<double,8> _c; //!< |r2|^2
		///@}

		
		/** @defgroup SplitFuncs Splitting function pointers
		 */
		///@{
		///@}

		/** @name LO splitting functions
		 *  @ingroup SplitFuncs
		 */
		///@{
		std::shared_ptr<SplittingFunction> _P0ns;
		std::shared_ptr<SplittingFunction> _P0qq;
		std::shared_ptr<SplittingFunction> _P0qg;
		std::shared_ptr<SplittingFunction> _P0gq;
		std::shared_ptr<SplittingFunction> _P0gg;
		///@}

		/** @name NLO splitting functions
		 *  @ingroup SplitFuncs
		 */
		///@{
		std::shared_ptr<SplittingFunction> _P1nsp;
		std::shared_ptr<SplittingFunction> _P1nsm;
		std::shared_ptr<SplittingFunction> _P1qq;
		std::shared_ptr<SplittingFunction> _P1qg;
		std::shared_ptr<SplittingFunction> _P1gq;
		std::shared_ptr<SplittingFunction> _P1gg;
		///@}

		/** @name NNLO splitting functions
		 *  @ingroup SplitFuncs
		 */
		///@{
		std::shared_ptr<SplittingFunction> _P2nsp;
		std::shared_ptr<SplittingFunction> _P2nsm;
		std::shared_ptr<SplittingFunction> _P2nsv;
		std::shared_ptr<SplittingFunction> _P2qq;
		std::shared_ptr<SplittingFunction> _P2qg;
		std::shared_ptr<SplittingFunction> _P2gq;
		std::shared_ptr<SplittingFunction> _P2gg;
		///@}

		/**
		 * @name NNLO splitting functions cache
		 * @ingroup SplitFuncs
		 */
		///@{
		std::vector<double> _P2nsp_cache;
		std::vector<double> _P2nsm_cache;
		std::vector<double> _P2nsv_cache;
		std::vector<double> _P2qq_cache;
		std::vector<double> _P2qg_cache;
		std::vector<double> _P2gq_cache;
		std::vector<double> _P2gg_cache;
		///@}

		/** @name NNNLO splitting functions
		 *  @ingroup SplitFuncs
		 */
		///@{
		std::shared_ptr<SplittingFunction> _P3nsp;
		std::shared_ptr<SplittingFunction> _P3nsm;
		std::shared_ptr<SplittingFunction> _P3nsv;
		///@}

		/**
		 * @name N3LO splitting functions cache
		 * @ingroup SplitFuncs
		 */
		///@{
		std::vector<double> _P3nsp_cache;
		std::vector<double> _P3nsm_cache;
		std::vector<double> _P3nsv_cache;
		///@}


		/** @name Operator Matrix Elements
		 */
		///@{
		std::shared_ptr<OpMatElem> _A2ns;
		std::shared_ptr<OpMatElem> _A2gq;
		std::shared_ptr<OpMatElem> _A2gg;
		std::shared_ptr<OpMatElem> _A2hq;
		std::shared_ptr<OpMatElem> _A2hg;
		///@}


		typedef std::chrono::high_resolution_clock Clock;
		


		/** @defgroup RecRels Recursion relations
		 */
		///@{
		///@}

		/** @name Singlet recursion relations
		 *  @ingroup RecRels
		 */
		///@{
		double RecRelS_1(std::vector<double> const& S,
						 uint k,
						 std::shared_ptr<SplittingFunction> P);
		double RecRelS_2(std::vector<double> const& S1,
						 std::vector<double> const& S2,
						 uint k,
						 std::shared_ptr<SplittingFunction> P0,
						 std::shared_ptr<SplittingFunction> P1);
		double RecRelS_3(std::vector<double> const& S1,
						 std::vector<double> const& S2,
						 std::vector<double> const& S3,
						 uint k,
						 std::shared_ptr<SplittingFunction> P0,
						 std::shared_ptr<SplittingFunction> P1,
						 std::shared_ptr<SplittingFunction> P2);
		///@}

		/** @name LO recursion relations
		 *  @ingroup RecRels
		 */
		///@{
		double RecRelLO(std::vector<double> const& A,
						uint k,
						std::shared_ptr<SplittingFunction> P);
		///@}

		/** @name NLO recursion relations
		 *  @ingroup RecRels
		 */
		///@{
		double RecRelNLO_1(std::vector<double> const& b,
						   uint k,
						   std::shared_ptr<SplittingFunction> P);
		double RecRelNLO_2(std::vector<double> const& B,
						   uint k,
						   std::shared_ptr<SplittingFunction> P);
		///@}

		/** @name NNLO recursion relations
		 *  @ingroup RecRels
		 */
		///@{
		double RecRelNNLO_1(std::vector<double> const& C,
							uint k,
							std::shared_ptr<SplittingFunction> P);
		double RecRelNNLO_2(std::vector<double> const& C,
							uint k,
							std::shared_ptr<SplittingFunction> P);
		double RecRelNNLO_3(std::vector<double> const& C,
							uint k,
							std::shared_ptr<SplittingFunction> P);
		///@}

		/** @name NNNLO recursion relations
		 *  @ingroup RecRels
		 */
		///@{
		double RecRelN3LO_1(std::vector<double> const& D,
							uint k,
							std::shared_ptr<SplittingFunction> P);
		double RecRelN3LO_2(std::vector<double> const& D,
							uint k,
							std::shared_ptr<SplittingFunction> P1,
							std::shared_ptr<SplittingFunction> P2,
							std::shared_ptr<SplittingFunction> P3);
		double RecRelN3LO_3(std::vector<double> const& D,
							uint k,
							std::shared_ptr<SplittingFunction> P1,
							std::shared_ptr<SplittingFunction> P2,
							std::shared_ptr<SplittingFunction> P3);
		double RecRelN3LO_4(std::vector<double> const& D,
							uint k,
							std::shared_ptr<SplittingFunction> P1,
							std::shared_ptr<SplittingFunction> P2,
							std::shared_ptr<SplittingFunction> P3);
		///@}

		
	public:
		/** @name Constructors/destructors
		 */
		///@{
		/** @brief Main constructor to initialize evolution
			@param order: perturbative order
		 *  @param grid: initialized grid object
		 *  @param Qf: final energy to evaluate to
		 *  @param initial_dist: initial distribution for masses, pdfs, and alpha0
		 */
	    DGLAPSolver(const uint order, Grid const& grid, const double Qf,
					std::shared_ptr<Distribution> initial_dist);
		~DGLAPSolver();
		///@}


		
		/** @brief Returns a const-ref to the \f$\alpha_s\f$ object.
		 */
		inline AlphaS const& GetAlphaS() const { return _alpha_s; }

		/** @brief Returns a const-ref to the @a Grid object.
		 */
		inline Grid const& GetGrid() const { return _grid; }


		/** @brief Main function to evolve coefficients.
		 */
		void Evolve();


		/** @brief Outputs a datafile in the new format
		 *  (very similar to old format but with cpp functions)
		 *  @param filepath: Name of output datafile
		*/
		void SetOutputDataFileNew(std::string const& filepath);

		/** @brief Outputs a datafile in the original Candia format
		 *  @param filepath: Name of output datafile
		*/
		void SetOutputDataFileOld(std::string const& filepath);
		

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

		/** @brief Computes (from the calculated coefficients)
		 */
		void SetupFinalDistributions();


		/** @brief evolves the singlet coefficients
		 */
		void EvolveSinglet();

		/** @brief evolves the non-singlet coefficients
		 */
		void EvolveNonSinglet();

		/** @brief After computing coefficients, do resummation to a tabulated energy value
		 *  @param Q: the tabulated energy value to evolve to
		 */
		void Resum(const double Q);

		/** @brief Finish resumming to the next threshold
		 */
		void ResumThreshold();

		/** Handles the matching conditions/transition operator matrix elements
		 */
		void HeavyFlavorTreatment();


		// Some debug-related stuff
	private:
		uint _output_file_index;

		// void OutputLOCoefficients();

		std::ofstream _debug_file;
	};

	

	
}


#endif // __CANDIA_HPP
