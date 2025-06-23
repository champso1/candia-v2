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
#include "Candia-v2/Math.hpp"

#include <memory>
#include <sstream>


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

		/** @name NNNLO splitting functions
		 *  @ingroup SplitFuncs
		 */
		///@{
		std::shared_ptr<SplittingFunction> _P3nsp;
		std::shared_ptr<SplittingFunction> _P3nsm;
		std::shared_ptr<SplittingFunction> _P3nss;
		///@}


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
		double RecRelNLO_2(std::vector<double> const& B1,
						   std::vector<double> const& B2,
						   uint k,
						   std::shared_ptr<SplittingFunction> P);
		///@}

		/** @name NNLO recursion relations
		 *  @ingroup RecRels
		 *  @note NOT IMPLEMENTED YET
		 */
		///@{
		double RecRelNNLO_1();
		double RecRelNNLO_2();
		double RecRelNNLO_3();
		///@}

		/** @name NNNLO recursion relations
		 *  @ingroup RecRels
		 *  @note NOT IMPLEMENTED YET
		 */
		///@{
		double RecRelN3LO_1();
		double RecRelN3LO_2();
		double RecRelN3LO_3();
		double RecRelN3LO_4();
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
		 *  (very similar but with cpp functions)
		 *  @param filepath: Name of output datafile
		*/
		void OutputDataFileNew(std::string const& filepath);

		/** @brief Outputs a datafile in the original Candia format
		 *  @param filepath: Name of output datafile
		*/
		void OutputDataFileOld(std::string const& filepath);
		

	private:
		/** @brief initializes coefficient arrays with distributions
		 */
		void SetInitialConditions();
		
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
		 */
		void Resum();

		/** @brief Finish resumming to the next threshold
		 */
		void ResumThreshold();
		
	};

	

	
}


#endif // __CANDIA_HPP
