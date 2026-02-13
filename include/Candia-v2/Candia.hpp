/**
 *  @file Candia.hpp
 *  @brief Contains the @a DGLAPSolver class which handles the actual evolution of the distributions
 */

#ifndef __CANDIA_HPP
#define __CANDIA_HPP

#include <concepts>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Grid.hpp"
#include "Candia-v2/AlphaS.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/ArrayGrid.hpp"
#include "Candia-v2/Options.hpp"

namespace Candia2
{
	/**
	 *  @brief groups all of the options/flags together into a single struct
	 */
	struct DGLAPOptions final
	{
		bool use_nnlo_matching_conditions_at_n3lo{true}; //!< switch for whether to use nnlo matching at n3lo (for benchmarking purposes)
		bool disable_heavy_flavor_matching{false}; //!< switch for whether to use matching at all for the heavy flavors
		bool use_n3lo_heavyquark_asymmetry{true}; //!< use new OME from arXiv:2512.13508
		bool use_truncated_nonsinglet_sol{false}; //!< whether to use the truncated ansatz in the non-singlet sector, as opposed to the exact solution
	};

	/**
	 *  @brief Performs the evolution of parton distributions.
	 */
	class DGLAPSolver : public OptionsBase<DGLAPOptions>
	{
	private:
		uint _order{}; //!< perturbative order
		Grid _grid; //!< main @a Grid object
		double _Qf{}; //!< final energy value to evolve to
		AlphaS _alpha_s; //!< main @a AlphaS object
		double _mur2_muf2{}; //!< \f$mu_r^2/mu_f^2\f$
		double _log_mur2_muf2{}; //!< \f$\log(\frac{\mu_r^2}{\mu_f^2})\f$
		double _log_muf2_mur2{}; //!< \f$\log(\frac{\mu_f^2}{\mu_r^2})\f$
		bool _is_scale_difference{}; //!< convenience flag for if \f$\mu_r \neq \mu_f\f$

		// TODO: remove the nfi and nff variables, they shouldn't be needed anymore
		// since they are stored inside the alpha_s object
		uint _nf{}; //!< current active number of massless flavors
		uint _nfi{}; //!< minimum number based on initial evolution and provided quark masses
		uint _nff{}; //!< final based on final evolution and provided quark masses
		double _alpha0{}; //!< initial alpha_s in an interval
		double _alpha1{}; //!< final alpha_s in an interval

		uint _iterations{}; //!< number of singlet/non-singlet iterations
		uint _trunc_idx{}; //!< number of additional singlet truncated iterations

		MultiDimArrayGrid_t<2> _A{}; //!< LO coeffs
		MultiDimArrayGrid_t<3> _B{}; //!< NLO coeffs
		MultiDimArrayGrid_t<4> _C{}; //!< NNLO coeffs
		MultiDimArrayGrid_t<5> _D{}; //!< N3LO coeffs
	    MultiDimArrayGrid_t<3> _S{}; //!< singlet coeffs
		MultiDimArrayGrid_t<3> _S_NS{}; //!< non-singlet coeffs (truncated)
		std::vector<ArrayGrid> _F{}; //!< final distributions

		std::array<double,8> _r1{}; //!< real solution to N3LO quadratic
		std::array<double,8> _b{};  //!< \f$-2*\mathrm{Re}[r_2]\f$
		std::array<double,8> _c{};  //!< \f$|r_2|^2\f$

		std::map<std::string_view, std::unique_ptr<Expression>> _expressions{}; //!< list of internal stores Expression objects
		/**
		 *  @brief Stores a unique pointer of the requested expression internally
		 *  @param name The name to associate the expression with
		 *  @param args Any arguments to pass to the constructor of the expression
		 */
		template <typename TExpr, typename... TExprArgs>
		requires (std::derived_from<TExpr, Expression>)
		void createExpression(std::string_view name, TExprArgs&&... args)
		{
			std::unique_ptr<Expression> ptr = std::make_unique<TExpr>(std::forward<TExprArgs>(args)...);
			_expressions.emplace(std::make_pair(name, std::move(ptr)));
		}

		/**
		 *  @brief returns the expression object with the given name
		 *  @param name The name associated with the expression object
		 */
		Expression& getExpression(std::string_view name);

		/** 
		 *  @brief Helper function to create/load all expression objects for the provided order.
		 */
		void loadAllExpressions();

	public:
		/**
		 *  @brief Constructs a DGLAPSolver object
		 *  @param order The perturbative order
		 *  @param grid reference to a @a Grid object
		 *  @param alpha_s (const) reference to an @a AlphaS object
		 *  @param Qf the final energy to evolve to
		 *  @param iterations Number of iterations to complete (the outer, s-index)
		 *  @param trunc_idx Number of terms in the singlet expansion to truncate at
		 *  @param initial_dist (const) reference to a @a Distribution object
		 *  @param mur2_muf2 The ratio of \f$\mu_r^2/\mu_f^2\f$
		 */
	    DGLAPSolver(
			uint order, Grid& grid, AlphaS const& alpha_s,
			double Qf, uint iterations, uint trunc_idx,
			Distribution const& initial_dist,
			double mur2_muf2 = 1.0);
		~DGLAPSolver(); //!< default destructor

		/** getter for the @a AlphaS object */
		inline AlphaS const& getAlphaS() const { return _alpha_s; }
		/** getter for the @a Grid object */
		inline Grid const& getGrid() const { return _grid; }

		/**
		 *  @brief Performs the full evolution.
		 */
		std::vector<ArrayGrid> const& evolve();

	private:
		/**
		 *  @brief Sets the initial values of the coefficients using the distribution
		 *  @param dist The initial distribution
		 */
		void setInitialConditions(Distribution const& dist);
		/**
		 *  @brief Sets the non-singlet/singlet distributions that are actually evolved rather than the raw quark distributions
		 */
		void setupCoefficients();
		/**
		 *  @brief Does the opposite of @a setInitialConditions, i.e. retrieves the raw quark dists from the special evolution ones
		 *  @param resum_ns the set of resummed non-singlet distributions
		 *  @param resum_singlet the set of resummed singlet distributions
		 *  @param resum[out] the final set of all (non-singlet and singlet) fixed distributions
		 */
		void fixDistributions(
			std::vector<ArrayGrid>& resum_ns, 
			std::vector<ArrayGrid>& resum_singlet,
			std::vector<ArrayGrid>& resum);

		/** @brief takes the default exact coefficients (A, B, ...) and sets up S to contain all necessary info */
		void setupTruncatedDistributions();

		void evolveSinglet(std::reference_wrapper<std::vector<ArrayGrid>> arr, double L1);
		/**
		 *  @brief Evolves the non-singlet distributions
		 *  @param arr Reference to the array in which to place the resummed results,
		 *  which will be different if resumming to the final energy vs a threshold one
		 *  @param L1 The LO logarithmic coefficient
		 *  @param L2 The NLO logarithmic coefficient
		 *  @param L3 The NNLO logarithmic coefficient (actually an arctan)
		 *  @param L4 The N3LO logarithmic coefficient
		 */
		void evolveNonSinglet(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			double L1, double L2, double L3, double L4);

		/**
		 *  @brief Evolves the non-singlet distributions with the truncated ansatz
		 *  @param arr Reference to the array in which to place the resummed results,
		 *  which will be different if resumming to the final energy vs a threshold one
		 *  @param L1 The LO logarithmic coefficient
		 */
		void evolveNonSingletTrunc(std::vector<ArrayGrid>& arr, double L1);
		

#if ENABLE_THREADING
		/**
		 *  @brief Evolves the singlet distributions
		 *  @param arr Reference to the array in which to place the resummed results,
		 *  which will be different if resumming to the final energy vs a threshold one
		 *  @param L1 the logarithmic coefficient
		 */
		void evolveSingletThreaded(
			std::reference_wrapper<std::vector<ArrayGrid>> arr,
			double L1);
		/**
		 *  @brief Evolves the non-singlet distributions
		 *  @param arr Reference to the array in which to place the resummed results,
		 *  which will be different if resumming to the final energy vs a threshold one
		 *  @param L1 The LO logarithmic coefficient
		 *  @param L2 The NLO logarithmic coefficient
		 *  @param L3 The NNLO logarithmic coefficient (actually an arctan)
		 *  @param L4 The N3LO logarithmic coefficient
		 */
		void evolveNonSingletThreaded(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			double L1, double L2, double L3, double L4);
#endif

		/**
		 *  @brief Performs the heavy-flavor matching for the quark/gluon distributions
		 */
		void heavyFlavorTreatment();

		/**
		 *  @defgroup hfthelpers Helpers for Heavy Flavor Treatment
		 *  @{
		 */
		/**
		 *  @brief relation 1 for NNLO HFT
		 *  @param c NNLO array grid (quark/anti-quark)
		 *  @param k current grid index
		 *  @param q quark array grid
		 */
		void HFT_NNLO1(ArrayGrid& c, uint k, ArrayGrid& q);
		/**
		 *  @brief relation 2 for NNLO HFT
		 *  @param g gluon arraygrid
		 *  @param qp q^(+) arraygrid
		 *  @param k current grid index
		 */
		void HFT_NNLO2(ArrayGrid& g, ArrayGrid& qp, uint k);
		/**
		 *  @brief relation 3 for NNLO HFT
		 *  @param g gluon arraygrid
		 *  @param qp q^(+) arraygrid
		 *  @param k current grid index
		 *  @param qh heavy quark array grid
		 *  @param qhb anti-heavy quark array grid
		 */
		void HFT_NNLO3(ArrayGrid& g, ArrayGrid& qp, uint k, ArrayGrid& qh, ArrayGrid& qhb);

		/**
		 *  @brief relation 1 for N3LO HFT
		 *  @param q quark arraygrid
		 *  @param qb anti-quark arraygrid
		 *  @param j distribution index
		 *  @param k current grid index
		 *  @param SP pre-calculated piece independent of @a j
		 *  @param qh the output q array to place the results. not the heavy quark dist, but a suitable enough name
		 */
		void HFT_N3LO1(ArrayGrid& q, ArrayGrid& qb, uint j, uint k, double SP, ArrayGrid& qh);
		/**
		 *  @brief relation 2 for N3LO HFT
		 *  @param q quark arraygrid
		 *  @param qb anti-quark arraygrid
		 *  @param j distribution index
		 *  @param k current grid index
		 *  @param SP pre-calculated piece independent of @a j
		 *  @param qhb the output qbar array to place the results. not the heavy quark dist, but a suitable enough name
		 */
		void HFT_N3LO2(ArrayGrid& q, ArrayGrid& qb, uint j, uint k, double SP, ArrayGrid& qhb);
		/**
		 *  @brief relation 3 for N3LO HFT
		 *  @param g gluon arraygrid
		 *  @param qp q^(+) arraygrid
		 *  @param k current grid index
		 */
		void HFT_N3LO3(ArrayGrid& g, ArrayGrid& qp, uint k);
		/**
		 *  @brief relation 4 for N3LO HFT
		 *  @param g gluon arraygrid
		 *  @param qp q^(+) arraygrid
		 *  @param qminus q^(-) arraygrid
		 *  @param k current grid index
		 *  @param qh the heavy quark array to place the results
		 *  @param qhb the heavy quark bar array to place the results
		 */
		void HFT_N3LO4(ArrayGrid& g, ArrayGrid& qp, ArrayGrid& qminus, uint k, ArrayGrid& qh, ArrayGrid& qhb);
		/** @} */


		/**
		 *  @defgroup singlethelpers Multi-Thread Singlet Helper Functions
		 *  @{
		 */
		/**
		 *  @brief Performs the NLO evolution of the coefficients
		 *  @param t truncation index
		 *  @param thread_idx a unique index passed to the function to ensure safe access to a gsl workspace
		 *  @param min minimum grid index
		 *  @param max maximum grid index
		 */
		void _mt_EvolveDistributions_S_NLO(uint t, int thread_idx, uint min, uint max);
		/**
		 *  @brief Performs the NNLO evolution of the coefficients
		 *  @param t truncation index
		 *  @param thread_idx a unique index passed to the function to ensure safe access to a gsl workspace
		 *  @param min minimum grid index
		 *  @param max maximum grid index
		 */
		void _mt_EvolveDistributions_S_NNLO(uint t, int thread_idx, uint min, uint max);
		/**
		 *  @brief Performs the N3LO evolution of the coefficients
		 *  @param t truncation index
		 *  @param thread_idx a unique index passed to the function to ensure safe access to a gsl workspace
		 *  @param min minimum grid index
		 *  @param max maximum grid index
		 */
		void _mt_EvolveDistributions_S_N3LO(uint t, int thread_idx, uint min, uint max);
		/** @} */

		/**
		 *  @defgroup nonsinglethelpers Multi-Thread Non-Singlet Helper Functions
		 *  @{
		 */
		/**
		 *  @brief LO non-singlet multi-threaded helper routine
		 *  @param arr reference to set of dists to place resummation results into
		 *  @param j distribution index
		 *  @param L1 log term
		 */
		void _mt_EvolveDistribution_NS_LO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, double L1);
		/**
		 *  @brief NLO non-singlet multi-threaded helper routine
		 *  @param arr reference to set of dists to place resummation results into
		 *  @param j distribution index
		 *  @param P1 the name of the NS piece of the P1 splitting function, e.g. the plus or minus piece
		 
		 *  @param L array of logarithm terms up to NLO
		 */
	    void _mt_EvolveDistribution_NS_NLO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::string const& P1, std::array<double, 2> const& L);
		/**
		 *  @brief NNLO non-singlet multi-threaded helper routine
		 *  @param arr reference to set of dists to place resummation results into
		 *  @param j distribution index
		 *  @param P array of names for the pieces of the NS splitting functions, e.g. the plus, minus, or valence pieces
		 *  @param L array of logarithm terms up to NNLO
		 */
	    void _mt_EvolveDistribution_NS_NNLO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::array<std::string, 2> const& P, std::array<double, 3> const& L);
		/**
		 *  @brief N3LO non-singlet multi-threaded helper routine
		 *  @param arr reference to set of dists to place resummation results into
		 *  @param j distribution index
		 *  @param P array of names for the pieces of the NS splitting functions, e.g. the plus, minus, or valence pieces
		 *  @param L array of logarithm terms up to N3LO
		 */
	    void _mt_EvolveDistribution_NS_N3LO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::array<std::string, 3> const& P, std::array<double, 4> const& L);
		/** @} */

		/** @defgroup recrels Recursion Relations */

		/**
		 *  @defgroup singletrecrels Singlet Recursion Relations
		 *  @ingroup recrels
		 *  @{
		 */
		/**
		 *  @brief LO singlet recursion relation: \f$S^i_{n+1}(x) = -\frac{2}{\beta_0}[P \otimes S^i_n](x)\f$
		 *  @param S arraygrid
		 *  @param k grid index
		 *  @param P LO splitting function
		 */
		double recrelS_1(
			ArrayGrid& S,
			uint k,
			Expression& P);
		/**
		 *  @brief NLO singlet recursion relation
		 *  @param S_i arraygrid corresponding to \f$S^i\f$
		 *  @param S_im1 arraygrid corresponding to \f$S^{i-1}\f$
		 *  @param k grid index
		 *  @param P0 LO splitting function
		 *  @param P1 NLO splitting function
		 */
		double recrelS_2(
			ArrayGrid& S_i, ArrayGrid& S_im1,
			uint k,
			Expression& P0, Expression& P1);
		/**
		 *  @brief NNLO singlet recursion relation
		 *  @param S_i arraygrid corresponding to \f$S^i\f$
		 *  @param S_im1 arraygrid corresponding to \f$S^{i-1}\f$
		 *  @param S_im2 arraygrid corresponding to \f$S^{i-2}\f$
		 *  @param k grid index
		 *  @param P0 LO splitting function
		 *  @param P1 NLO splitting function
		 *  @param P2 NNLO splitting function
		 */
		double recrelS_3(
			ArrayGrid& S_i, ArrayGrid& S_im1, ArrayGrid& S_im2,
			uint k,
			Expression& P0, Expression& P1, Expression& P2);
		/**
		 *  @brief N3LO singlet recursion relation
		 *  @param S_i arraygrid corresponding to \f$S^i\f$
		 *  @param S_im1 arraygrid corresponding to \f$S^{i-1}\f$
		 *  @param S_im2 arraygrid corresponding to \f$S^{i-2}\f$
		 *  @param S_im3 arraygrid corresponding to \f$S^{i-3}\f$
		 *  @param k grid index
		 *  @param P0 LO splitting function
		 *  @param P1 NLO splitting function
		 *  @param P2 NNLO splitting function
		 *  @param P3 N3LO splitting function
		 */
		double recrelS_4(
			ArrayGrid& S_i, ArrayGrid& S_im1, ArrayGrid& S_im2, ArrayGrid& S_im3,
			uint k,
			Expression& P0, Expression& P1, Expression& P2, Expression& P3);
		/** @} */

		/**
		 *  @defgroup lononsingletrecrels LO Non-Singlet Recursion Relations
		 *  @ingroup recrels
		 *  @{
		 */
		/**
		 *  @brief performs the LO recursion relation \f$A_{n+1}(x) = -\frac{2}{\beta_0}[P \otimes A_n](x)\f$
		 *  @param A the LO arraygrid
		 *  @param k the grid index
		 *  @param P0 the LO splitting function
		 */
		double recrelLO(
			ArrayGrid& A,
			uint k,
			Expression& P0);
		/** @} */

		/**
		 *  @defgroup nlononsingletrecrels NLO Non-Singlet Recursion Relations
		 *  @ingroup recrels
		 *  @{
		 */
		/**
		 *  @brief performs the NLO recursion relation equivalent to the LO one
		 *  @param B the NLO arraygrid
		 *  @param k the grid index
		 *  @param P0 the LO splitting function
		 */
		double recrelNLO_1(
			ArrayGrid& B,
			uint k,
			Expression& P0);
		/**
		 *  @brief performs the 2nd NLO recursion relation
		 *  @param B the NLO arraygrid
		 *  @param k the grid index
		 *  @param P0 the LO splitting function
		 *  @param P1 the NLO splitting function
		 */
		double recrelNLO_2(
			ArrayGrid& B,
			uint k,
			Expression& P0, Expression& P1);
		/** @} */

		/**
		 *  @defgroup nnlononsingletrecrels NNLO Non-Singlet Recursion Relations
		 *  @ingroup recrels
		 *  @{
		 */
		/**
		 *  @brief performs the NNLO recursion relation equivalent to the LO one
		 *  @param C the NNLO arraygrid
		 *  @param k the grid index
		 *  @param P0 the LO splitting function
		 */
		double recrelNNLO_1(
			ArrayGrid& C,
			uint k,
			Expression& P0);
		/**
		 *  @brief performs the 2nd NNLO recursion relation
		 *  @param C the NNLO arraygrid
		 *  @param k the grid index
		 *  @param P0 the LO splitting function
		 *  @param P1 the NLO splitting function
		 *  @param P2 the NNLO splitting function
		 */
		double recrelNNLO_2(
			ArrayGrid& C,
			uint k,
			Expression& P0, Expression& P1, Expression& P2);
		/**
		 *  @brief performs the 3rd NNLO recursion relation
		 *  @param C the NNLO arraygrid
		 *  @param k the grid index
		 *  @param P0 the LO splitting function
		 *  @param P1 the NLO splitting function
		 */
		double recrelNNLO_3(
			ArrayGrid& C,
			uint k,
			Expression& P0, Expression& P1);
		/** @} */

		/**
		 *  @defgroup n3lononsingletrecrels N3LO Non-Singlet Recursion Relations
		 *  @ingroup recrels
		 *  @{
		 */
		/**
		 *  @brief performs the N3LO recursion relation equivalent to the LO one
		 *  @param D the NNLO arraygrid
		 *  @param k the grid index
		 *  @param P0 the LO splitting function
		 */
		double recrelN3LO_1(
			ArrayGrid& D,
			uint k,
			Expression& P0);
		/**
		 *  @brief performs the 2nd N3LO recursion relation
		 *  @param D the NNLO arraygrid
		 *  @param k the grid index
		 *  @param P0 the LO splitting function
		 *  @param P1 the NLO splitting function
		 *  @param P2 the NNLO splitting function
		 *  @param P3 the N3LO splitting function
		 */
		double recrelN3LO_2(
			ArrayGrid& D,
			uint k,
			Expression& P0, Expression& P1, Expression& P2, Expression& P3);
		/**
		 *  @brief performs the 3nd N3LO recursion relation
		 *  @param D the NNLO arraygrid
		 *  @param k the grid index
		 *  @param P0 the LO splitting function
		 *  @param P1 the NLO splitting function
		 *  @param P2 the NNLO splitting function
		 *  @param P3 the N3LO splitting function
		 */
		double recrelN3LO_3(
			ArrayGrid& D,
			uint k,
			Expression& P0, Expression& P1, Expression& P2, Expression& P3);
		/**
		 *  @brief performs the 4th N3LO recursion relation
		 *  @param D the NNLO arraygrid
		 *  @param k the grid index
		 *  @param P0 the LO splitting function
		 *  @param P1 the NLO splitting function
		 *  @param P2 the NNLO splitting function
		 *  @param P3 the N3LO splitting function
		 */
		double recrelN3LO_4(
			ArrayGrid& D,
			uint k,
			Expression& P0, Expression& P1, Expression& P2, Expression& P3);
		/** @} */
	};
}


#endif // __CANDIA_HPP
