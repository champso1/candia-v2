#ifndef __CANDIA_HPP
#define __CANDIA_HPP


#include "Candia-v2/Common.hpp"
#include "Candia-v2/Grid.hpp"
#include "Candia-v2/AlphaS.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/FuncArrayGrid.hpp"

#include <concepts>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

namespace Candia2
{	
	class DGLAPSolver
	{
	private:
		uint _order{}; //!< perturbative order
		Grid _grid; //!< main @a Grid object
		double _Qf{}; //!< final energy value to evolve to
		AlphaS _alpha_s; //!< main @a AlphaS object
		double _mur2_muf2{}; //!< mu_r^2/mu_f^2
		double _log_mur2_muf2{}; //!< log(_mur2_muf2)

		// TODO: remove the nfi and nff variables, they shouldn't be needed anymore
		// since they are stored inside the alpha_s object
		uint _nf{}; //!< current active number of massless flavors
		uint _nfi{}; //!< minimum number based on initial evolution and provided quark masses
		uint _nff{}; //!< final based on final evolution and provided quark masses
		double _alpha0{}; //!< initial alpha_s in a threshold
		double _alpha1{}; //!< final alpha_s in a threshold

		uint _iterations{}; //!< number of singlet/non-singlet iterations
		uint _trunc_idx{}; //!< number of additional singlet truncated iterations

		MultiDimArrayGrid_t<2> _A2{}; //!< LO coeffs
		MultiDimArrayGrid_t<3> _B2{}; //!< NLO coeffs
		MultiDimArrayGrid_t<4> _C2{}; //!< NNLO coeffs
		MultiDimArrayGrid_t<5> _D2{}; //!< N3LO coeffs
	    MultiDimArrayGrid_t<3> _S2{}; //!< singlet coeffs
		std::vector<ArrayGrid> _F2{}; //!< final distributions

		std::array<double,8> _r1{}; //!< one real solution to N3LO quadratic
		std::array<double,8> _b{};  //!< -2*Re[r2]
		std::array<double,8> _c{};  //!< |r2|^2

		// TODO: make this a CMAKE definition or some macro of some kind
		// that way we can avoid linking with pthread or the Windows equivalent
		// which would be nice for debug builds, since at the moment
		// I believe we link with it no matter what
		bool _multi_thread;

		// for purposes of comparing with benchmarks,
		// this flag lets one enable whether or not to use
		// the n3lo matching conditions in the n3lo evolution
		bool _use_n3lo_matching_conditions;

		std::map<std::string_view, std::unique_ptr<Expression>> _expressions{};
		template <typename TExpr, typename... TExprArgs>
		requires (std::derived_from<TExpr, Expression>)
		void createExpression(std::string_view name, Grid const& grid, TExprArgs&&... args)
		{
			(void)grid; // TODO: remove this
			std::unique_ptr<Expression> ptr = std::make_unique<TExpr>(std::forward<TExprArgs>(args)...);
			_expressions.emplace(std::make_pair(name, std::move(ptr)));
		}
		Expression& getExpression(std::string_view name);
		void loadAllExpressions();

	public:
	    DGLAPSolver(
			uint order, Grid const& grid, AlphaS const& alpha_s,
			double Qf, uint iterations, uint trunc_idx,
			Distribution const& initial_dist,
			double kr = 1.0, bool multi_thread = false);
		~DGLAPSolver();
		
		inline AlphaS const& getAlphaS() const { return _alpha_s; }
		inline Grid const& getGrid() const { return _grid; }
		inline void useNNLOMatchingAtN3LO() { _use_n3lo_matching_conditions = false;};
		
		void setEvolutionVariables(uint iterations, uint trunc_idx);

		auto evolve() -> decltype(_F2);

	private:
		void setInitialConditions(Distribution const& dist);
		void fillSplittingFunctionCaches();
		void setupCoefficients();
		void fixDistributions(
			bool resum_tab, bool resum_threshold, 
			std::vector<ArrayGrid>& temp_arr, 
			std::vector<ArrayGrid>& temp_arr_singlet);

		void evolveSinglet(std::reference_wrapper<std::vector<ArrayGrid>> arr, double L1);
		void evolveNonSinglet(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			double L1, double L2, double L3, double L4);
		void evolveNonSingletThreaded(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			double L1, double L2, double L3, double L4);

	
		void heavyFlavorTreatment();
		void HFT_NNLO1(ArrayGrid& c, uint k, ArrayGrid& q);
		void HFT_NNLO2(ArrayGrid& g, ArrayGrid& qp, uint k);
		void HFT_NNLO3(ArrayGrid& g, ArrayGrid& qp, uint k, ArrayGrid& qh, ArrayGrid& qhb);
		
		void HFT_N3LO1(ArrayGrid& q, ArrayGrid& qb, uint j, uint k, double SP);
		void HFT_N3LO2(ArrayGrid& q, ArrayGrid& qb, uint j, uint k, double SP);
		void HFT_N3LO3(ArrayGrid& g, ArrayGrid& qp, uint k);
		void HFT_N3LO4(ArrayGrid& g, ArrayGrid& qp, uint k);

		void _mt_EvolveDistribution_NS_LO  (
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, double L1);
	    void _mt_EvolveDistribution_NS_NLO (
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::string const& P1, std::array<double, 2> const& L);
	    void _mt_EvolveDistribution_NS_NNLO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::array<std::string, 2> const& P, std::array<double, 3> const& L);
	    void _mt_EvolveDistribution_NS_N3LO(
			std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::array<std::string, 3> const& P, std::array<double, 4> const& L);


		double recrelS_1(
			ArrayGrid& S,
			uint k,
			Expression& P);
		double recrelS_2(
			ArrayGrid& S1, ArrayGrid& S2,
			uint k,
			Expression& P0, Expression& P1);
		double recrelS_3(
			ArrayGrid& S1, ArrayGrid& S2, ArrayGrid& S3,
			uint k,
			Expression& P0, Expression& P1, Expression& P2);
		double recrelS_4(
			ArrayGrid& S1, ArrayGrid& S2, ArrayGrid& S3, ArrayGrid& S4,
			uint k,
			Expression& P0, Expression& P1, Expression& P2, Expression& P3);

		double recrelLO(
			ArrayGrid& A,
			uint k,
			Expression& P0);
		double recrelNLO_1(
			ArrayGrid& B,
			uint k,
			Expression& P0);
		double recrelNLO_2(
			ArrayGrid& B,
			uint k,
			Expression& P0, Expression& P1);

		double recrelNNLO_1(
			ArrayGrid& C,
			uint k,
			Expression& P0);
		double recrelNNLO_2(
			ArrayGrid& C,
			uint k,
			Expression& P0, Expression& P1, Expression& P2);
		double recrelNNLO_3(
			ArrayGrid& C,
			uint k,
			Expression& P0, Expression& P1);

		double recrelN3LO_1(
			ArrayGrid& D,
			uint k,
			Expression& P0);
		double recrelN3LO_2(
			ArrayGrid& D,
			uint k,
			Expression& P0, Expression& P1, Expression& P2, Expression& P3);
		double recrelN3LO_3(
			ArrayGrid& D,
			uint k,
			Expression& P0, Expression& P1, Expression& P2, Expression& P3);
		double recrelN3LO_4(
			ArrayGrid& D,
			uint k,
			Expression& P0, Expression& P1, Expression& P2, Expression& P3);
	};
}


#endif // __CANDIA_HPP
