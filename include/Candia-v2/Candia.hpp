/**
 *  @file Candia.hpp
 *  @brief Contains the @a DGLAPSolver class which handles the actual evolution of the distributions
 */

#pragma once

#include <concepts>
#include <memory>
#include <utility>
#include <vector>
#include <optional>

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Couplings.hpp"
#include "Candia-v2/Expression.hpp"
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
		bool use_nnlo_matching_conditions_at_n3lo{false}; //!< switch for whether to use nnlo matching at n3lo (for benchmarking purposes)
		bool disable_heavy_flavor_matching{false}; //!< switch for whether to use matching at all for the heavy flavors
		bool use_n3lo_heavyquark_asymmetry{true}; //!< use new OME from arXiv:2512.13508		 
		bool use_fortran_nnlo_splitfuncs{false}; //!< whether to use the fortran versions for the nnlo splitting functions or the C++-translated ones
		bool use_fortran_n3lo_splitfuncs{false}; //!< whether to use the fortran versions for the n3lo splitting functions or the C++-translated ones
		bool try_qed{false}; //!< whether to try and include QED effects
	};

	/**
	 *  @brief Performs the evolution of parton distributions.
	 */
	class DGLAPSolver : public OptionsBase<DGLAPOptions>
	{
	public:
		using options_type = DGLAPOptions;
	private:
		uint _order{}; //!< perturbative order
		Grid _grid; //!< main @a Grid object
		double _Qf{}; //!< final energy value to evolve to
		AlphaS _alpha_s; //!< main @a AlphaS object
		std::optional<AlphaQED> _alpha_qed; //!< main @a AlphaQED object
		double _mur2_muf2{}; //!< \f$mu_r^2/mu_f^2\f$
		double _log_mur2_muf2{}; //!< \f$\log(\frac{\mu_r^2}{\mu_f^2})\f$
		double _log_muf2_mur2{}; //!< \f$\log(\frac{\mu_f^2}{\mu_r^2})\f$
		bool _is_scale_difference{}; //!< convenience flag for if \f$\mu_r \neq \mu_f\f$

		// TODO: remove the nfi and nff variables, they shouldn't be needed anymore
		// since they are stored inside the alpha_s object
		uint _nf{}; //!< current active number of massless flavors
		const uint _nl{3}; //!< number of active leptons (set to 3 by default for now)
		uint _nfi{}; //!< minimum number based on initial evolution and provided quark masses
		uint _nff{}; //!< final based on final evolution and provided quark masses
		double _alpha0{}; //!< initial alpha_s in an interval
		double _alpha1{}; //!< final alpha_s in an interval
		double _alphaqed0{}; //!< initial alpha_qed in an interval
		double _alphaqed1{}; //!< final alpha_qed in an interval

		inline std::vector<uint> getEvolutionIndices()
		{
			std::vector<uint> idxs{};
			switch (_order) {
				case 0:
				case 1: {
					for (uint j=13; j<=12+_nf; ++j)
						idxs.emplace_back(j);
					for (uint j=32; j<=30+_nf; ++j)
						idxs.emplace_back(j);
				}; break;
				case 2:
				case 3: {
					for (uint j=26; j<=24+_nf; ++j)
						idxs.emplace_back(j);
					for (uint j=32; j<=30+_nf; ++j)
						idxs.emplace_back(j);
					idxs.emplace_back(25);
				} break;
				default: throw std::runtime_error("unreachable");
			}
			return idxs;
		}

		Distribution const& _initial_dist; //!< reference to initial distribution

		uint _iterations{}; //!< number of singlet/non-singlet iterations
		uint _trunc_idx{}; //!< number of additional singlet truncated iterations

	    ArrayGridN<3> _A{}; //!< LO coeffs
		ArrayGridN<4> _B{}; //!< NLO coeffs
		ArrayGridN<5> _C{}; //!< NNLO coeffs
		ArrayGridN<6> _D{}; //!< N3LO coeffs
		ArrayGridN<4> _S{}; //!< singlet coeffs
		ArrayGridN<4> _S_NS{}; //!< non-singlet coeffs (truncated)
		std::vector<ArrayGrid> _F{}; //!< final distributions

		ArrayGridN<4> _A_QED{}; //!< LO QCD/LO QED non-singlet coeffs
		ArrayGridN<4> _S_QED{}; //!< LO QCD/LO QED singlet coeffs

		std::array<double,8> _r1{}; //!< real solution to N3LO quadratic
		std::array<double,8> _b{};  //!< \f$-2*\mathrm{Re}[r_2]\f$
		std::array<double,8> _c{};  //!< \f$|r_2|^2\f$

		std::array<std::unique_ptr<Expression>, static_cast<uint>(ExprName::Count)> _expressions{};
		template <typename TExpr, typename... TExprArgs>
		requires (std::derived_from<TExpr, Expression>)
		inline void createExpression(ExprName name, TExprArgs&&... args)
		{
		    _expressions[static_cast<uint>(name)] = std::make_unique<TExpr>(std::forward<TExprArgs>(args)...);
		}
	public:
		/**
		 *  @brief returns the expression object with the given name
		 *  @param name The name associated with the expression object
		 */
		inline Expression& getExpression(ExprName name)
		{
			return *_expressions[static_cast<uint>(name)];
		}

	private:
		enum class EvolType
		{
			None,
			Exact,
			Truncated
		};
		EvolType _evol_type{EvolType::None};

		/** @defgroup coeffgetters Helpers for  */
		/**
		 *  @ingroup coeffgetters
		 *  @{
		 */
		inline double& getSingletCoeffValue(uint j, uint k) {
			assert(_evol_type != EvolType::None, "invalid evolution type (found None)");
			if (_evol_type == EvolType::Exact)
				return _S(0,j,0,k);
			else
				return _S_NS(0,j,0,k);
		}
		inline ArrayGridView getSingletCoeffArray(uint j) {
			assert(_evol_type != EvolType::None, "invalid evolution type (found None)");
			if (_evol_type == EvolType::Exact)
				return _S(0,j,0);
			else
				return _S_NS(0,j,0);
		}
		
	    inline double& getNonSingletCoeffValue(uint j, uint k) {
			assert(_evol_type != EvolType::None, "invalid evolution type (found None)");
			if (_evol_type == EvolType::Truncated)
				return _S_NS(0,j,0,k);
			switch (_order) {
				case 0: return _A(j,0,k); break;
				case 1: return _B(j,0,0,k); break;
				case 2: return _C(j,0,0,0,k); break;
				case 3: return _D(j,0,0,0,0,k); break;
				default: throw std::runtime_error("unreachable");
			}
		}

		inline ArrayGridView getNonSingletCoeffArray(uint j) {
			assert(_evol_type != EvolType::None, "invalid evolution type (found None)");
			if (_evol_type == EvolType::Truncated)
				return _S_NS(0,j,0);
			switch (_order) {
				case 0: return _A(j,0); break;
				case 1: return _B(j,0,0); break;
				case 2: return _C(j,0,0,0); break;
				case 3: return _D(j,0,0,0,0); break;
				default: throw std::runtime_error("unreachable");
			}
		}
		
		/** 
		 *  @brief Helper function to create/load all expression objects for the provided order.
		 */
		void loadAllExpressions();


	private:
		std::array<P3ApproxType, static_cast<uint>(ExprName::Count)> _p3_approx_types{
			([](){
				std::array<P3ApproxType, static_cast<uint>(ExprName::Count)> p3approxs{};
				std::ranges::fill(p3approxs, P3ApproxType::ImodAvg);
				p3approxs[static_cast<uint>(ExprName::P3nsm)] = P3ApproxType::Exact;
				p3approxs[static_cast<uint>(ExprName::P3nsp)] = P3ApproxType::Exact;
				p3approxs[static_cast<uint>(ExprName::P3nsv)] = P3ApproxType::Exact;
				return p3approxs;
			})()
		};

	public:
		inline void setP3ApproximationTypes(std::vector<std::pair<ExprName,P3ApproxType>> types)
		{
			for (auto const& [first, second] : types)
				_p3_approx_types[static_cast<uint>(first)] = second;
		}


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
			uint order,
			Grid& grid,
			AlphaS const& alpha_s,
			double Qf, uint iterations, uint trunc_idx,
			Distribution const& initial_dist,
			double mur2_muf2 = 1.0);
		
		/**
		 *  @brief Constructs a DGLAPSolver object with QED effects
		 *  @param order The perturbative order
		 *  @param grid reference to a @a Grid object
		 *  @param alpha_s (const) reference to an @a AlphaS object
		 *  @param alpha_qed (const) reference to an @a AlphaQED object
		 *  @param Qf the final energy to evolve to
		 *  @param iterations Number of iterations to complete (the outer, s-index)
		 *  @param trunc_idx Number of terms in the singlet expansion to truncate at
		 *  @param initial_dist (const) reference to a @a Distribution object
		 *  @param mur2_muf2 The ratio of \f$\mu_r^2/\mu_f^2\f$
		 */
	    DGLAPSolver(
			uint order,
			Grid& grid,
			AlphaS const& alpha_s, AlphaQED const& alpha_qed,
			double Qf, uint iterations, uint trunc_idx,
			Distribution const& initial_dist,
			double mur2_muf2 = 1.0);
		~DGLAPSolver(); //!< default destructor

		/** getter for the @a AlphaS object */
		inline AlphaS const& getAlphaS() const { return _alpha_s; }
		/** getter for the @a Grid object */
		inline Grid& getGrid() { return _grid; }

		/**
		 *  @brief Performs the full evolution using the exact ansatz in the NS sector
		 */
		inline std::vector<ArrayGrid> const& evolve()
		{
			log(LOG_DEBUG, "DGLAPSolver", "Using exact ansatz for NS sector");
			return _evolve_function(EvolType::Exact);
		}

		/**
		 *  @brief Performs the full evolution using the truncated ansatz in the NS sector
		 */
		inline std::vector<ArrayGrid> const& evolveTrunc()
		{
			log(LOG_DEBUG, "DGLAPSolver", "Using exact ansatz for NS sector");
			return _evolve_function(EvolType::Truncated);
		}

		/**
		 *  @brief returns a vector of some subtraction PDFs as in EQ.27, Eq.38 in arXiv:2410.03876 [hep-ph]
		 */
		std::vector<ArrayGrid> calculateSubtractionPDFs();
	private:
		std::vector<ArrayGrid> const& _evolve_function(EvolType evoltype);
	private:
		/**
		 *  @brief Sets the non-singlet/singlet distributions that are actually evolved rather than the raw quark distributions
		 */
		void setupCoefficients();
		/**
		 *  @brief Sets up the QED basis
		 */
		void setupQEDCoefficients();

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

		/**
		 *  @brief Does the opposite of @a setInitialConditions, i.e. retrieves the raw quark dists from the special evolution ones (QED)
		 *  @param resum_ns the set of resummed non-singlet distributions
		 *  @param resum_singlet the set of resummed singlet distributions
		 *  @param resum[out] the final set of all (non-singlet and singlet) fixed distributions
		 */
		void fixDistributionsQED(
			std::vector<ArrayGrid>& resum_ns, 
			std::vector<ArrayGrid>& resum_singlet,
			std::vector<ArrayGrid>& resum);

	    /**
		 * same as @a fixDistributions(), but doesn't zero parts of the arrays in the event no evolution has occurred yet
		 * this is because the resum_{singlet,ns} arrays won't have been set up
		 */
		void fixDistributionsForce(std::vector<ArrayGrid>& resum);

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
		void evolveNonSinglet(std::reference_wrapper<std::vector<ArrayGrid>> arr, double L1, double L2, double L3, double L4);

		/**
		 *  @brief Evolves the non-singlet distributions with the truncated ansatz
		 *  @param arr Reference to the array in which to place the resummed results,
		 *  which will be different if resumming to the final energy vs a threshold one
		 *  @param L1 The LO logarithmic coefficient
		 */
		void evolveNonSingletTrunc(std::reference_wrapper<std::vector<ArrayGrid>> arr, double L1);

		/**
		 *  @brief Performs the PDF evolution including QED effects
		 *  @param arr Reference to the array in which to place the resummed results,
		 *  which will be different if resumming to the final energy vs a threshold one
		 *  @param L1 The LO logarithmic coefficient
		 */
		void evolveQED(
			std::reference_wrapper<std::vector<ArrayGrid>> arr_singlet,
			std::reference_wrapper<std::vector<ArrayGrid>> arr_ns,
			double L0QED, double L0QCD);
		
		/**
		 *  @brief Performs the heavy-flavor matching for the quark/gluon distributions
		 */
		void heavyFlavorTreatment();

		/**
		 *  @defgroup hfthelpers Helpers for Heavy Flavor Treatment
		 *  @{
		 */
		void HFT_NNLO1(ArrayGridView c, uint k, ArrayGridView q);
		void HFT_NNLO2(ArrayGridView g, ArrayGridView qp, uint k);
		void HFT_NNLO3(ArrayGridView g, ArrayGridView qp, uint k, ArrayGridView qh, ArrayGridView qhb);

		void HFT_N3LO1(ArrayGridView q, ArrayGridView qb, uint k, double SP, ArrayGridView qh);
		void HFT_N3LO2(ArrayGridView q, ArrayGridView qb, uint k, double SP, ArrayGridView qhb);
		void HFT_N3LO3(ArrayGridView g, ArrayGridView qp, uint k);
		void HFT_N3LO4(ArrayGridView g, ArrayGridView qp, ArrayGridView qminus, uint k, ArrayGridView qh, ArrayGridView qhb);
		/** @} */


		/**
		 *  @defgroup nonsinglethelpers Multi-Thread Non-Singlet Helper Functions
		 *  @{
		 */
		void _mt_EvolveDistribution_NS_LO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, double L1);
	    void _mt_EvolveDistribution_NS_NLO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, ExprName P1, std::array<double, 2> const& L);
	    void _mt_EvolveDistribution_NS_NNLO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::array<ExprName, 2> const& P, std::array<double, 3> const& L);

	    void _mt_EvolveDistribution_NS_N3LO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::array<ExprName, 3> const& P, std::array<double, 4> const& L);
		/** @} */

		/**
		 *  @defgroup nonsinglettrunchelpers Multi-Thread Non-Singlet Helper Functions (Truncated Ansatz)
		 *  @{
		 */
		void _mt_EvolveDistribution_NST_LO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr,
			uint j, double L1);
	    void _mt_EvolveDistribution_NST_NLO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr,
			uint j, ExprName p1, double L1);
	    void _mt_EvolveDistribution_NST_NNLO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr,
			uint j, ExprName p1, ExprName p2, double L1);
	    void _mt_EvolveDistribution_NST_N3LO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr,
			uint j, ExprName p1, ExprName p2, ExprName p3, double L1);
		/** @} */

		/** @defgroup recrels Recursion Relations */
		/** @defgroup recrelhelpers Recursion Relation Helpers */

		/**
		 *  @ingroup recrelhelpers
		 *  @{
		 */
		double shift_p1(Expression& P1, Expression& P0, ArrayGridView A, uint k);
		double shift_p2(Expression& P2, Expression& P1, Expression& P0, ArrayGridView A, uint k);
		double shift_p3(Expression& P3, Expression& P2, Expression& P1, Expression& P0, ArrayGridView A, uint k);
		/** @} */
		
		/**
		 *  @defgroup singletrecrels Singlet Recursion Relations
		 *  @ingroup recrels
		 *  @{
		 */
		double recrelS_1(ArrayGridView S, uint k, Expression& P);
		double recrelS_2(
			ArrayGridView S_i, ArrayGridView S_im1,
			uint k,
			Expression& P0, Expression& P1);
		double recrelS_3(
			ArrayGridView S_i, ArrayGridView S_im1, ArrayGridView S_im2,
			uint k,
			Expression& P0, Expression& P1, Expression& P2);
		double recrelS_4(
			ArrayGridView S_i, ArrayGridView S_im1, ArrayGridView S_im2, ArrayGridView S_im3,
			uint k,
			Expression& P0, Expression& P1, Expression& P2, Expression& P3);
		/** @} */

		/**
		 *  @defgroup lononsingletrecrels LO Non-Singlet Recursion Relations
		 *  @ingroup recrels
		 *  @{
		 */
		double recrelLO(ArrayGridView A, uint k, Expression& P0);
		/** @} */

		/**
		 *  @defgroup nlononsingletrecrels NLO Non-Singlet Recursion Relations
		 *  @ingroup recrels
		 *  @{
		 */
		double recrelNLO_1(
			ArrayGridView B,
			uint k,
			Expression& P0);
		double recrelNLO_2(
			ArrayGridView B,
			uint k,
			Expression& P0, Expression& P1);
		/** @} */

		/**
		 *  @defgroup nnlononsingletrecrels NNLO Non-Singlet Recursion Relations
		 *  @ingroup recrels
		 *  @{
		 */
		double recrelNNLO_1(
			ArrayGridView C,
			uint k,
			Expression& P0);
		double recrelNNLO_2(
			ArrayGridView C,
			uint k,
			Expression& P0, Expression& P1, Expression& P2);
		double recrelNNLO_3(
			ArrayGridView C,
			uint k,
			Expression& P0, Expression& P1);
		/** @} */

		/**
		 *  @defgroup n3lononsingletrecrels N3LO Non-Singlet Recursion Relations
		 *  @ingroup recrels
		 *  @{
		 */
		double recrelN3LO_1(
			ArrayGridView D,
			uint k,
			Expression& P0);
		double recrelN3LO_2(
			ArrayGridView D,
			uint k,
			Expression& P0, Expression& P1, Expression& P2, Expression& P3);
		double recrelN3LO_3(
			ArrayGridView D,
			uint k,
			Expression& P0, Expression& P1, Expression& P2, Expression& P3);
		double recrelN3LO_4(
			ArrayGridView D,
			uint k,
			Expression& P0, Expression& P1, Expression& P2, Expression& P3);
		/** @} */
	};
}
