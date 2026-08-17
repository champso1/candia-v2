// Candia.cpp

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/Expression.hpp"
#include "Candia-v2/Grid.hpp"
#include "Candia-v2/SplittingFn.hpp"
#include "Candia-v2/SplittingFnQED.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"
#include "Candia-v2/Math.hpp"
#include <algorithm>


// PDF indices
// 
// 0      gluon
// 1-6    quarks
// 7-12   antiquarks
// 13-18  q_i^-
// 19-24  q_i^+
// 25     q^(-)
// 26-30  q_{NS,1i}^(-)
// 31     q^(+)
// 32-36  q_{NS,1i}^(+)


namespace
{
	inline constexpr char const* CANDIA_OPENING_TEXT = 
		"============================================================\n"
    	"| \033[1mCandia-v2\033[0m — DGLAP evolution up to \033[1mN³LO\033[0m\n"
    	"| Based on Candia (C version)\n"
    	"| \033[2mPlease cite \033[4marXiv:2512.22667\033[0m\033[2m for this version and\033[0m\n"
		"| \033[2m\033[4marXiv:0803.0462\033[0m\033[2m for Candia-v1\033[0m\n"
    	"| \033[2m(Hampson, Guzzi)\033[0m\n"
		"============================================================\n";
		inline constexpr char const* CANDIA_CLOSING_TEXT = 
		"============================================================\n"
    	"| \033[1mCandia-v2\033[0m — DGLAP evolution up to \033[1mN³LO\033[0m\n"
    	"| Thanks for using this program!\n"
        "| \033[2mPlease cite \033[4marXiv:2512.22667\033[0m\033[2m for this version and\033[0m\n"
		"| \033[2m\033[4marXiv:0803.0462\033[0m\033[2m for Candia-v1\033[0m\n"
    	"| \033[2m(Hampson, Guzzi)\033[0m\n"
		"============================================================\n";

}

namespace Candia2
{

	DGLAPSolver::DGLAPSolver(
		uint order,
		Grid& grid,
		AlphaS const& alpha_s,
		double Qf, uint iterations, uint trunc_idx,
		Distribution const& initial_dist,
		double mur2_muf2) 
		: _order{order},  _grid{grid}, _Qf{Qf},
		  _alpha_s{alpha_s},
		  _mur2_muf2{mur2_muf2}, _log_mur2_muf2{std::log(mur2_muf2)}, _log_muf2_mur2{-std::log(mur2_muf2)},
		  _is_scale_difference{mur2_muf2 != 1.0},
		  _initial_dist{initial_dist},
		  _iterations{iterations}, _trunc_idx{trunc_idx},
		  _A({DISTS, 2, grid.size()}),
		  _B({DISTS, 2, iterations, grid.size()}),
		  _C({DISTS, 2, iterations, iterations, grid.size()}),
		  _D({DISTS, 2, iterations, iterations, iterations, grid.size()}),
		  _S({trunc_idx+1, 2, 2, grid.size()}),
		  _S_NS{},
		  _F(DISTS, ArrayGrid(grid.size()))
	{
		log(::CANDIA_OPENING_TEXT);

		log(LOG_DEBUG, "DGLAPSolver::DGLAPSolver()", "Evolving with log(mu_R / mu_F) = log({:.1}) = {:.4}.", _mur2_muf2, _log_mur2_muf2);

		if (_order == 0 && _trunc_idx != 0) {
			_trunc_idx = 0; // LO has exact singlet solution, do not add additional terms
			log(LOG_WARNING, "DGLAP", "Specified value of the truncation index ({}) will be set to zero.", _trunc_idx);
		}

		_r1[1] = -0.965105642503553;
		_b[1] = -2.0 * 0.1629296392275606;
		_c[1] = std::pow(0.1629296392275606, 2) + std::pow(0.9535744823175397, 2);

		_r1[2] = -1.0315080774348302;
		_b[2] = -2.0 * 0.18523659836580222;
		_c[2] = std::pow(0.18523659836580222, 2) + std::pow(1.0299109343730084, 2);
				
		_r1[3] = -1.1120253073038324;
		_b[3] = -2.0 * 0.2214224000789979;
		_c[3] = std::pow(0.2214224000789979, 2) + std::pow(1.131077812338495, 2);

		_r1[4] = -1.2090185772488318;
		_b[4] = -2.0 * 0.2867586032664649;
		_c[4] = std::pow(0.2867586032664649, 2) + std::pow(1.272794345339416, 2);
				
		_r1[5] = -1.3205899823870375;
		_b[5] = -2.0 * 0.42477034063852415;
		_c[5] = std::pow(0.42477034063852415, 2) + std::pow(1.4854822725151384, 2);
				
		_r1[6] = -1.4277979273114205;
		_b[6] = -2.0 * 0.7964970177083996;
		_c[6] = std::pow(0.7964970177083996, 2) + std::pow(1.816809978388145, 2);

		auto func = [](std::array<double, 8> const& a) -> std::string {
			auto view =
				std::views::iota(1) | std::views::take(6)
				| std::views::transform([&a](int i){ return a[i]; });
			return vec_to_str(view);
		};
		log(LOG_DEBUG, "DGLAPSolver::DGLAPSolver()", "r1 array: {}", func(_r1));
		log(LOG_DEBUG, "DGLAPSolver::DGLAPSolver()", "b  array: {}", func(_b));
		log(LOG_DEBUG, "DGLAPSolver::DGLAPSolver()", "c  array: {}", func(_c));
	}

	DGLAPSolver::DGLAPSolver(
		uint order,
		Grid& grid,
		AlphaS const& alpha_s, AlphaQED const& alpha_qed,
		double Qf, uint iterations, uint trunc_idx,
		Distribution const& initial_dist,
		double mur2_muf2) 
		: _order{order},  _grid{grid}, _Qf{Qf},
		  _alpha_s{alpha_s}, _alpha_qed{alpha_qed},
		  _mur2_muf2{mur2_muf2}, _log_mur2_muf2{std::log(mur2_muf2)}, _log_muf2_mur2{-std::log(mur2_muf2)},
		  _is_scale_difference{mur2_muf2 != 1.0},
		  _initial_dist{initial_dist},
		  _iterations{iterations}, _trunc_idx{trunc_idx}
		  // _A({DISTS, 2, grid.size()}),
		  // _B({DISTS, 2, iterations, grid.size()}),
		  // _C({DISTS, 2, iterations, iterations, grid.size()}),
		  // _D({DISTS, 2, iterations, iterations, iterations, grid.size()}),
		  // _S({trunc_idx+1, 2, 2, grid.size()}),
		  // _S_NS{},
		  // _F(DISTS, ArrayGrid(grid.size()))
	{
		log(::CANDIA_OPENING_TEXT);

		log(LOG_DEBUG, "DGLAPSolver::DGLAPSolver()", "Evolving with log(mu_R / mu_F) = log({:.1}) = {:.4}.", _mur2_muf2, _log_mur2_muf2);

		if (_order == 0 && _trunc_idx != 0) {
			_trunc_idx = 0; // LO has exact singlet solution, do not add additional terms
			log(LOG_WARNING, "DGLAP", "Specified value of the truncation index ({}) will be set to zero.", _trunc_idx);
		}

		_r1[1] = -0.965105642503553;
		_b[1] = -2.0 * 0.1629296392275606;
		_c[1] = std::pow(0.1629296392275606, 2) + std::pow(0.9535744823175397, 2);

		_r1[2] = -1.0315080774348302;
		_b[2] = -2.0 * 0.18523659836580222;
		_c[2] = std::pow(0.18523659836580222, 2) + std::pow(1.0299109343730084, 2);
				
		_r1[3] = -1.1120253073038324;
		_b[3] = -2.0 * 0.2214224000789979;
		_c[3] = std::pow(0.2214224000789979, 2) + std::pow(1.131077812338495, 2);

		_r1[4] = -1.2090185772488318;
		_b[4] = -2.0 * 0.2867586032664649;
		_c[4] = std::pow(0.2867586032664649, 2) + std::pow(1.272794345339416, 2);
				
		_r1[5] = -1.3205899823870375;
		_b[5] = -2.0 * 0.42477034063852415;
		_c[5] = std::pow(0.42477034063852415, 2) + std::pow(1.4854822725151384, 2);
				
		_r1[6] = -1.4277979273114205;
		_b[6] = -2.0 * 0.7964970177083996;
		_c[6] = std::pow(0.7964970177083996, 2) + std::pow(1.816809978388145, 2);

		auto func = [](std::array<double, 8> const& a) -> std::string {
			auto view =
				std::views::iota(1) | std::views::take(6)
				| std::views::transform([&a](int i){ return a[i]; });
			return vec_to_str(view);
		};
		log(LOG_DEBUG, "DGLAPSolver::DGLAPSolver()", "r1 array: {}", func(_r1));
		log(LOG_DEBUG, "DGLAPSolver::DGLAPSolver()", "b  array: {}", func(_b));
		log(LOG_DEBUG, "DGLAPSolver::DGLAPSolver()", "c  array: {}", func(_c));
	}

	DGLAPSolver::~DGLAPSolver()
	{
		log(LOG_INFO, "DGLAP", "Exiting...");
		log(::CANDIA_CLOSING_TEXT);
	}

	void DGLAPSolver::loadAllExpressions()
    {
        createExpression<P0ns>(ExprName::P0ns);
        createExpression<P0qq>(ExprName::P0qq);
        createExpression<P0qg>(ExprName::P0qg);
        createExpression<P0gq>(ExprName::P0gq);
        createExpression<P0gg>(ExprName::P0gg);

		createExpression<P1nsm>(ExprName::P1nsm);
		createExpression<P1nsp>(ExprName::P1nsp);
		createExpression<P1qq>(ExprName::P1qq);
		createExpression<P1qg>(ExprName::P1qg);
		createExpression<P1gq>(ExprName::P1gq);
		createExpression<P1gg>(ExprName::P1gg);

		if (options.use_fortran_nnlo_splitfuncs) {
			log(LOG_DEBUG, "DGLAPSolver::loadAllExpressions()", "Loading Fortran versions of P2 splitting functions");
			createExpression<mvv_p2::P2nsm>(ExprName::P2nsm);
			createExpression<mvv_p2::P2nsp>(ExprName::P2nsp);
			createExpression<mvv_p2::P2nsv>(ExprName::P2nsv);
			createExpression<mvv_p2::P2qq>(ExprName::P2qq);
			createExpression<mvv_p2::P2qg>(ExprName::P2qg);
			createExpression<mvv_p2::P2gq>(ExprName::P2gq);
			createExpression<mvv_p2::P2gg>(ExprName::P2gg);
		} else {
			log(LOG_DEBUG, "DGLAPSolver::loadAllExpressions()", "Loading C++ versions of P2 splitting functions");
			createExpression<P2nsm>(ExprName::P2nsm);
			createExpression<P2nsp>(ExprName::P2nsp);
			createExpression<P2nsv>(ExprName::P2nsv);
			createExpression<P2qq>(ExprName::P2qq);
			createExpression<P2qg>(ExprName::P2qg);
			createExpression<P2gq>(ExprName::P2gq);
			createExpression<P2gg>(ExprName::P2gg);
		}
			
		createExpression<A2ns>(ExprName::A2ns);
		createExpression<A2gq>(ExprName::A2gq);
		createExpression<A2gg>(ExprName::A2gg);
		createExpression<A2hq>(ExprName::A2hq);
		createExpression<A2hg>(ExprName::A2hg);
		auto it = std::ranges::find(_p3_approx_types, P3ApproxType::Exact);
		bool using_exact_p3ns = it != _p3_approx_types.end();
		
		if (options.use_fortran_n3lo_splitfuncs && using_exact_p3ns)
			log(LOG_ERROR, "DGLAPSolver::loadAllExpressions()", "Cannot use both the Fortran P3s and the exact (the Fortran versions are approximate)");
		
		if (options.use_fortran_n3lo_splitfuncs) {
			log(LOG_DEBUG, "DGLAPSolver::loadAllExpression()", "Loading Fortran versions of P3 splitting functions");
			createExpression<mvv_p3::P3nsm>(ExprName::P3nsm, _p3_approx_types[static_cast<uint>(ExprName::P3nsm)]);
			createExpression<mvv_p3::P3nsp>(ExprName::P3nsp, _p3_approx_types[static_cast<uint>(ExprName::P3nsp)]);
			createExpression<mvv_p3::P3nsv>(ExprName::P3nsv, _p3_approx_types[static_cast<uint>(ExprName::P3nsv)]);
			createExpression<mvv_p3::P3qq>(ExprName::P3qq, _p3_approx_types[static_cast<uint>(ExprName::P3qq)]);
			createExpression<mvv_p3::P3qg>(ExprName::P3qg, _p3_approx_types[static_cast<uint>(ExprName::P3qg)]);
			createExpression<mvv_p3::P3gq>(ExprName::P3gq, _p3_approx_types[static_cast<uint>(ExprName::P3gq)]);
			createExpression<mvv_p3::P3gg>(ExprName::P3gg, _p3_approx_types[static_cast<uint>(ExprName::P3gg)]);
		} else {
			createExpression<P3qq>(ExprName::P3qq, _p3_approx_types[static_cast<uint>(ExprName::P3qq)]);
			createExpression<P3qg>(ExprName::P3qg, _p3_approx_types[static_cast<uint>(ExprName::P3qg)]);
			createExpression<P3gq>(ExprName::P3gq, _p3_approx_types[static_cast<uint>(ExprName::P3gq)]);
			createExpression<P3gg>(ExprName::P3gg, _p3_approx_types[static_cast<uint>(ExprName::P3gg)]);

			if (_p3_approx_types[static_cast<uint>(ExprName::P3nsm)] == P3ApproxType::Exact) {
				log(LOG_DEBUG, "DGLAPSolver::loadAllExpressions()", "using exact P3nsm");
				createExpression<p3_exact::P3nsm>(ExprName::P3nsm);
			} else {
				log(LOG_DEBUG, "DGLAPSolver::loadAllExpressions()", "using approx P3nsm");
				createExpression<P3nsm>(ExprName::P3nsm, _p3_approx_types[static_cast<uint>(ExprName::P3nsm)]);
			}
			if (_p3_approx_types[static_cast<uint>(ExprName::P3nsp)] == P3ApproxType::Exact) {
				log(LOG_DEBUG, "DGLAPSolver::loadAllExpressions()", "using exact P3nsp");
				createExpression<p3_exact::P3nsp>(ExprName::P3nsp);
			} else {
				log(LOG_DEBUG, "DGLAPSolver::loadAllExpressions()", "using approx P3nsp");
				createExpression<P3nsp>(ExprName::P3nsp, _p3_approx_types[static_cast<uint>(ExprName::P3nsp)]);
			}
			if (_p3_approx_types[static_cast<uint>(ExprName::P3nsv)] == P3ApproxType::Exact) {
				log(LOG_DEBUG, "DGLAPSolver::loadAllExpressions()", "using exact P3nsv");
				createExpression<p3_exact::P3nsv>(ExprName::P3nsv);
			} else {
				log(LOG_DEBUG, "DGLAPSolver::loadAllExpressions()", "using approx P3nsv");
				createExpression<P3nsv>(ExprName::P3nsv, _p3_approx_types[static_cast<uint>(ExprName::P3nsv)]);
			}
		}

		createExpression<OpMatElemN3LO>(ExprName::A3nsm, ome::AqqQNSEven);
		createExpression<OpMatElemN3LO>(ExprName::A3nsp, ome::AqqQNSOdd);
		createExpression<OpMatElemN3LO>(ExprName::A3gq, ome::AgqQ);
		createExpression<OpMatElemN3LO>(ExprName::A3gg, ome::AggQ);
		createExpression<OpMatElemN3LO>(ExprName::A3hq, ome::AQqPS);
		createExpression<OpMatElemN3LO>(ExprName::A3hg, ome::AQg);
		createExpression<OpMatElemN3LO>(ExprName::A3psqq, ome::AqqQPS);
		createExpression<OpMatElemN3LO>(ExprName::A3sqg, ome::AqgQ);
		createExpression<OpMatElemN3LO>(ExprName::A3PSshq, ome::AQqPSs);

		log(LOG_DEBUG, "DGLAPSolver::loadAllExpressions()", "Using the following P3 approximation types:");
		for (uint i=static_cast<uint>(ExprName::P3nsm); i<=static_cast<uint>(ExprName::P3gg); ++i) {
			ExprName exprname = static_cast<ExprName>(i);
			P3ApproxType approxtype = _p3_approx_types[i];

			std::string exprname_str{}, approxtype_str{};
			switch (exprname) {
				case ExprName::P3nsm: exprname_str = "P3nsm"; break;
				case ExprName::P3nsp: exprname_str = "P3nsp"; break;
				case ExprName::P3nsv: exprname_str = "P3nsv"; break;
				case ExprName::P3qq: exprname_str = "P3qq"; break;
				case ExprName::P3qg: exprname_str = "P3qg"; break;
				case ExprName::P3gq: exprname_str = "P3gq"; break;
				case ExprName::P3gg: exprname_str = "P3gg"; break;
				default: throw std::runtime_error("unreachable");
			}
			switch (approxtype) {
				case P3ApproxType::Imod1: approxtype_str = "Imod1"; break;
				case P3ApproxType::Imod2: approxtype_str = "Imod2"; break;
				case P3ApproxType::ImodAvg: approxtype_str = "ImodAvg"; break;
				case P3ApproxType::Exact: approxtype_str = "Exact"; break;
			}
			log(LOG_DEBUG, "DGLAPSolver::loadAllExpressions()", "  - {} => {}", exprname_str, approxtype_str);
		}


		if (getOptions().try_qed) {
			log(LOG_DEBUG, "DGLAPSolver::loadAllExpressions()", "trying to load QED expressions");
			createExpression<P0ff>(ExprName::P0ff);
			createExpression<P0fy>(ExprName::P0fy);
			createExpression<P0yf>(ExprName::P0yf);
			createExpression<P0yy>(ExprName::P0yy);
		}
    }

    void DGLAPSolver::setupCoefficients()
    {
	    for (uint k=0; k<_grid.size(); k++) {
			for (uint j=13; j<=18; j++)
				getNonSingletCoeffValue(j, k) = getNonSingletCoeffValue(j-12, k)-getNonSingletCoeffValue(j-6, k);
		
			getNonSingletCoeffValue(25, k)=0.;
			for (uint j=13; j<=18; j++)
				getNonSingletCoeffValue(25, k) += getNonSingletCoeffValue(j, k);

			for (uint j=26; j<=30; j++)
				getNonSingletCoeffValue(j, k) = getNonSingletCoeffValue(13, k)-getNonSingletCoeffValue(j-12, k);

			for (uint j=19; j<=24; j++)
				getNonSingletCoeffValue(j, k) = getNonSingletCoeffValue(j-18, k)+getNonSingletCoeffValue(j-12, k);

			_S(0,1,0,k) = 0.0;
			for (uint j=19; j<=24; j++)
				_S(0,1,0,k) += getNonSingletCoeffValue(j, k);
		
			for (uint j=32; j<=36; j++)
				getNonSingletCoeffValue(j, k)=getNonSingletCoeffValue(19, k)-getNonSingletCoeffValue(j-12, k);
		}
    }

	void DGLAPSolver::setupQEDCoefficients()
    {
	    log(LOG_WARNING, "DGLAPSolver::setupQEDCoefficients()", "not implemented (nor needed, atm?)");
    }

	void DGLAPSolver::fixDistributions(
		std::vector<ArrayGrid>& resum_ns,
		std::vector<ArrayGrid>& resum_singlet,
		std::vector<ArrayGrid>& resum)
    {
        for (uint t=1; t<=_trunc_idx; ++t) {
			for (uint j=0; j<=1; ++j) {
				for (uint n=0; n<=1; ++n) 
					std::ranges::fill(_S(t,j,n), 0.0);
			}
        }

		if (_evol_type == EvolType::Truncated) {
			for (uint t=1; t<=_trunc_idx; ++t) {
				for (uint j=0; j<_F.size(); ++j) {
					for (uint n=0; n<=1; ++n) 
						std::ranges::fill(_S_NS(t,j,n), 0.0);
				}
			}
		}

		for (uint j=0; j<=1; ++j)
			std::ranges::copy(resum_singlet[j*31], resum[j*31].begin());
		switch (_order) {
			case 0:
			case 1: {
				for (uint j=13; j<=12+_nf; j++)
					std::ranges::copy(resum_ns[j], resum[j].begin());
				for (uint j=32; j<=30+_nf; j++)
					std::ranges::copy(resum_ns[j], resum[j].begin());
			} break;
			case 2:
			case 3: {
				for (uint j=26; j<=24+_nf; ++j)
					std::ranges::copy(resum_ns[j], resum[j].begin());
				for (uint j=32; j<=30+_nf; ++j)
					std::ranges::copy(resum_ns[j], resum[j].begin());
				std::ranges::copy(resum_ns[25], resum[25].begin());
			} break;
		}

        double Nf = static_cast<double>(_nf);
		for (uint k=0; k<_grid.size()-1;k++) {
			if (_order>=2) {
				resum[13](k) = resum[25](k);
				for (uint j=26; j<=24+_nf; j++)
					resum[13](k) += resum[j](k);
				resum[13](k) /= Nf;
				for (uint j=14; j<=12+_nf; j++)
					resum[j](k) = resum[13](k) - resum[j+12](k);
			}

			resum[19](k) = resum[31](k);
			for (uint j=32; j<=30+_nf; j++)
				resum[19](k) += resum[j](k);
			resum[19](k) /= Nf;

			for (uint j=20; j<=18+_nf; j++)
				resum[j](k) = resum[19](k) - resum[j+12](k);

			for (uint j=1; j<=_nf; j++) {
				resum[j](k)   =0.5*(resum[j+18](k) + resum[j+12](k));
				resum[j+6](k) =0.5*(resum[j+18](k) - resum[j+12](k));
			}

			if (_order<2) {
				resum[25](k)=0.0;
				for (uint j=13; j<=12+_nf; j++)
					resum[25](k) += resum[j](k);

				for (uint j=26; j<=24+_nf; j++)
					resum[j](k) = resum[13](k) - resum[j-12](k);
			}
		}
    }

	void DGLAPSolver::fixDistributionsQED(
		std::vector<ArrayGrid>& resum_ns,
		std::vector<ArrayGrid>& resum_singlet,
		std::vector<ArrayGrid>& resum)
    {
		auto num_dists = static_cast<uint>(QEDPartonIndices::COUNT);
		
	    auto sigmaud = static_cast<uint>(QEDPartonIndices::SIGMAUD);
		auto sigma = static_cast<uint>(QEDPartonIndices::SIGMA);
		auto gluon = static_cast<uint>(QEDPartonIndices::G);
		auto photon = static_cast<uint>(QEDPartonIndices::PHOTON);
		auto sigmal = static_cast<uint>(QEDPartonIndices::SIGMAL);
		std::array s_dists{sigmaud, sigma, gluon, photon, sigmal};

		auto ns_dists =
			std::views::iota(uint{0}, num_dists)
			| std::views::filter([&](uint j){ return std::ranges::find(s_dists, j) == s_dists.end(); });

		for (uint j : s_dists)
			std::ranges::copy(resum_singlet[j], resum[j].begin());
		for (uint j : ns_dists)
			std::ranges::copy(resum_ns[j], resum[j].begin());
    }

	void DGLAPSolver::fixDistributionsForce(std::vector<ArrayGrid>& resum)
	{
		for (uint j=0; j<=1; ++j)
			std::ranges::copy(getSingletCoeffArray(j), resum[j*31].begin());
		switch (_order) {
			case 0:
			case 1: {
				for (uint j=13; j<=12+_nf; j++)
					std::ranges::copy(getNonSingletCoeffArray(j), resum[j].begin());
				for (uint j=32; j<=30+_nf; j++)
					std::ranges::copy(getNonSingletCoeffArray(j), resum[j].begin());
			} break;
			case 2:
			case 3: {
				for (uint j=26; j<=24+_nf; ++j)
					std::ranges::copy(getNonSingletCoeffArray(j), resum[j].begin());
				for (uint j=32; j<=30+_nf; ++j)
					std::ranges::copy(getNonSingletCoeffArray(j), resum[j].begin());
				std::ranges::copy(getNonSingletCoeffArray(25), resum[25].begin());
			} break;
		}

        double Nf = static_cast<double>(_nf);
		for (uint k=0; k<_grid.size()-1;k++) {
			if (_order>=2) {
				resum[13](k)=resum[25](k);
				for (uint j=26; j<=24+_nf; j++)
					resum[13](k) += resum[j](k);
				resum[13](k) /= Nf;
				for (uint j=14; j<=12+_nf; j++)
					resum[j](k) = resum[13](k) - resum[j+12](k);
			}

			resum[19](k) = resum[31](k);
			for (uint j=32; j<=30+_nf; j++)
				resum[19](k) += resum[j](k);
			resum[19](k) /= Nf;

			for (uint j=20; j<=18+_nf; j++)
				resum[j](k) = resum[19](k) - resum[j+12](k);

			for (uint j=1; j<=_nf; j++) {
				resum[j](k)   = 0.5*(resum[j+18](k) + resum[j+12](k));
				resum[j+6](k) = 0.5*(resum[j+18](k) - resum[j+12](k));
			}

			if (_order<2) {
				resum[25](k)=0.0;
				for (uint j=13; j<=12+_nf; j++)
					resum[25](k) += resum[j](k);

				for (uint j=26; j<=24+_nf; j++)
					resum[j](k) = resum[13](k) - resum[j-12](k);
			}
		}
	}

	void DGLAPSolver::setupTruncatedDistributions()
	{
		log(LOG_DEBUG, "Grid", "Using truncated ansatz for non-singlet sector.");

		_S_NS.resize({_trunc_idx+1, DISTS, 2, _grid.size()});

		for (uint j=0; j<DISTS; ++j) {
			switch (_order) {
				case 0: std::ranges::copy(_A(j,0), _S_NS(0,j,0).begin()); break;
				case 1: std::ranges::copy(_B(j,0,0), _S_NS(0,j,0).begin()); break;
				case 2: std::ranges::copy(_C(j,0,0,0), _S_NS(0,j,0).begin()); break;
				case 3: std::ranges::copy(_D(j,0,0,0,0), _S_NS(0,j,0).begin()); break;
			}
		}
		
		_A.clear();
		_B.clear();
		_C.clear();
		_D.clear();
	}


	std::vector<ArrayGrid> const& DGLAPSolver::_evolve_function(EvolType evol_type)
	{
		_evol_type = evol_type;
		if (getOptions().try_qed && !_alpha_qed.has_value())
			log(LOG_ERROR, "DGLAPSolver::_evolve_function()", "wanting to do QED without having initialized alpha_qed");

		uint num_dists;
		if (getOptions().try_qed) {
			num_dists = static_cast<uint>(QEDPartonIndices::COUNT);
			_A_QED.resize({num_dists, 2, _iterations, _grid.size()});
			_S_QED.resize({num_dists, 2, _iterations, _grid.size()});

			auto s_accessor =
				[&](uint j, uint k) -> double& {
					return _S_QED(j,0,0,k);
				};
			auto ns_accessor =
				[&](uint j, uint k) -> double& {
					switch (_order) {
						case 0: return _A_QED(j,0,0,k); break;
						default: throw std::runtime_error("implement >LO for QED");
					}
				};
			_initial_dist.fillCoeffs(s_accessor, ns_accessor, _grid.points());
		} else {
			num_dists = static_cast<uint>(StandardPartonIndices::COUNT);
			auto s_accessor =
				[&](uint j, uint k) -> double& {
					return _S(0,j,0,k);
				};
			if (_evol_type == EvolType::Truncated) {
			    _S_NS.resize({_trunc_idx+1, num_dists, 2, _grid.size()});
				auto ns_accessor =
					[&](uint j, uint k) -> double& {
					    return _S_NS(0,j,0,k);
					};
				_initial_dist.fillCoeffs(s_accessor, ns_accessor, _grid.points());
			} else {
				switch (_order) {
					case 0: _A.resize({num_dists, 2, _grid.size()}); break;
					case 1: _B.resize({num_dists, 2, _iterations, _grid.size()}); break;
					case 2: _C.resize({num_dists, 2, _iterations, _iterations, _grid.size()}); break;
					case 3: _D.resize({num_dists, 2, _iterations, _iterations, _iterations, _grid.size()}); break;
					default: throw std::runtime_error("unreachable");
				}
				auto ns_accessor =
					[&](uint j, uint k) -> double& {
						switch (_order) {
							case 0: return _A(j,0,k); break;
							case 1: return _B(j,0,0,k); break;
							case 2: return _C(j,0,0,0,k); break;
							case 3: return _D(j,0,0,0,0,k); break;
							default: throw std::runtime_error("unreachable");
						}
					};
				_initial_dist.fillCoeffs(s_accessor, ns_accessor, _grid.points());
			}   
		}
	    
		
	    log(LOG_INFO, "DGLAP", "Evolving to {} flavors.", _alpha_s.nff());
		loadAllExpressions();
		
		std::vector<ArrayGrid> final_dists(num_dists, ArrayGrid(_grid.size()));

		// since we now only store two iterations at once,
		// we create these temporary arrays that will store the results of the resummation
		// since they were originally stored in the s=0 part that no longer exists
		// they are then copied to the _F final dists or the s=0 of the next iteration
		std::vector<ArrayGrid> resum_ns(num_dists, ArrayGrid(_grid.size()));
		std::vector<ArrayGrid> resum_singlet(num_dists, ArrayGrid(_grid.size()));
		std::vector<ArrayGrid> resum(num_dists, ArrayGrid(_grid.size()));
		
		bool performed_evolution = false;
		for (_nf=_alpha_s.nfi(); _nf<=_alpha_s.nff(); _nf++) {
			log(LOG_DEBUG, "DGLAP", "Setting nf={}", _nf);
			bool last_loop = _nf == _alpha_s.nff();

			log(LOG_DEBUG, "DGLAP", "Setting up distributions for evolution.");
			if (getOptions().try_qed)
				setupQEDCoefficients();
			else
				setupCoefficients();
			log(LOG_DEBUG, "DGLAP", "Finished setting up distributions for evolution.");

			// if the next mass is zero, we are already done
			// but this shouldn't really be hit
			// TODO: handle lepton masses?
			if (_alpha_s.masses(_nf+1) == 0.0) {
				log(LOG_WARNING, "DGLAP", "Next mass is zero. Quitting...");
				break;
			}

			// update all values
			_alpha_s.update(_nf);
			log(LOG_DEBUG, "DGLAP", "Loading relevant splitting function / OME values into cache");
			SplittingFunction::update(_nf, _alpha_s.beta0(), _log_muf2_mur2);
			OpMatElem::update(_log_muf2_mur2, _nf);
			for (auto& expr : _expressions)
				expr->preCalc();
			
			log(LOG_DEBUG, "DGLAP", "Retrieving values of alpha_s, and calculating all logarithm factors");
			bool resum_tab = _alpha_s.resumTabulated();
			bool resum_threshold = !resum_tab;
			_alpha0 = _alpha_s.post(_nf);
			_alpha1 = resum_tab ? 
				_alpha_s.evaluate(_alpha_s.masses(_nf), _Qf, _alpha0) :
				_alpha_s.pre(_nf+1);
			
			double beta0 = _alpha_s.beta0();
			double beta1 = _alpha_s.beta1();
			double beta2 = _alpha_s.beta2();
			double beta0qed{};
			double r1 = _r1[_nf];
			double b = _b[_nf];
			double c = _c[_nf];
			double L1 = std::log(_alpha1/_alpha0);
			double L1QED{};
			double L2{}, L3{}, L4{};
			if (_order == 1) {
				L2 = std::log((_alpha1*beta1 + 4.0*PI*beta0)
								/(_alpha0*beta1 + 4.0*PI*beta0));
			} else if (_order == 2) {
				L2 = std::log((16.0*PI_2*beta0 + 4.0*PI*_alpha1*beta1 + _alpha1*_alpha1*beta2)
							 /(16.0*PI_2*beta0 + 4.0*PI*_alpha0*beta1 + _alpha0*_alpha0*beta2));
				
				// analytic continuation for arctan
				double aux=4.0*beta0*beta2 - beta1*beta1;
				if (aux > 0) {
					L3 = std::atan(2.0*PI*(_alpha1-_alpha0)*std::sqrt(aux)
						          /(2.*PI*(8.*PI*beta0+(_alpha1+_alpha0)*beta1)+_alpha1*_alpha0*beta2)
					)/std::sqrt(aux);
				} else {
					log(LOG_WARNING, "DGLAP", "Encountered argument of square root that is <0. Doing analytic continuation.");
					L3 = std::atanh(2.*PI*(_alpha1-_alpha0)*std::sqrt(-aux)
						           /(2.*PI*(8.*PI*beta0+(_alpha1+_alpha0)*beta1)+_alpha1*_alpha0*beta2)
					)/std::sqrt(-aux);
				}
			} else if (_order == 3) {
				// NOTE: equation 2 and recursion relation 2 are determined from the log with the quadratic terms,
				// which were initially called L3, so I swapped them
				L3 = std::log((_alpha1 - r1)/(_alpha0 - r1));
				L2 = std::log((_alpha1*_alpha1 + b*_alpha1 + c)/(_alpha0*_alpha0 + b*_alpha0 + c));
				double aux = std::sqrt(-b*b + 4*c); // never negative, no need for analytic continuation
				L4 = std::atan((_alpha1-_alpha0)*aux / (2.0*_alpha0*_alpha1 + (_alpha0+_alpha1)*b + 2.0*c))/aux;
			}

			log(LOG_DEBUG, "DGLAP::evolve()", "Values of QCD log coeffs:");
			std::vector<std::pair<double, std::string_view>> coeffs{{L1, "L1"}, {L2, "L2"}, {L3, "L3"}, {L4, "L4"}};
			for (auto [x, xname] : coeffs)
				log(LOG_DEBUG, "DGLAP::evolve()", "  - {} = {: }", xname, x);
			
			log(LOG_DEBUG, "DGLAP", "Doing {} resummation", (resum_tab ? "tabulated" : "threshold" ));
			log(LOG_DEBUG, "DGLAP", "AlphaS: {} --> {}", _alpha0, _alpha1);
			
			if (getOptions().try_qed) {
				_alpha_qed.value().update(_nf, _nl);
				SplitFuncQED::update(_nf, _alpha_s.beta0(), _log_muf2_mur2, _nl);
				std::tie(_alphaqed0, _alphaqed1) = _alpha_qed.value().initFinalAlpha();
				L1QED = std::log(_alphaqed1/_alphaqed0);
				beta0qed = _alpha_qed.value().beta0();
				log(LOG_DEBUG, "DGLAP::evolve()", "Values of QED log coeffs:");
				log(LOG_DEBUG, "DGLAP::evolve()", "  - L1 = {: }", L1QED);
				log(LOG_DEBUG, "DGLAP::evolve()", "  - beta0 = {: }", beta0qed);
				log(LOG_DEBUG, "DGLAP", "AlphaQED: {} --> {}", _alphaqed0, _alphaqed1);
			}

			// only do evolution if alphas are different
			// (i.e. energy scales are different)
			if (_alpha0 != _alpha1) {
				performed_evolution = true;
				if (getOptions().try_qed) {
					log(LOG_DEBUG, "DGLAP", "Performing the evolution with QED effects");

					// singlet
					{
						auto& p0ns = getExpression(ExprName::P0ns);
						auto& p0qq = getExpression(ExprName::P0qq);
						auto& p0qg = getExpression(ExprName::P0qg);
						auto& p0gq = getExpression(ExprName::P0gq);
						auto& p0gg = getExpression(ExprName::P0gg);
						auto& p0ff = getExpression(ExprName::P0ff);
						auto& p0fy = getExpression(ExprName::P0fy);
						auto& p0yf = getExpression(ExprName::P0yf);
						auto& p0yy = getExpression(ExprName::P0yy);

						auto sigmaud = static_cast<uint>(QEDPartonIndices::SIGMAUD);
						auto sigma = static_cast<uint>(QEDPartonIndices::SIGMA);
						auto gluon = static_cast<uint>(QEDPartonIndices::G);
						auto photon = static_cast<uint>(QEDPartonIndices::PHOTON);
						auto sigmal = static_cast<uint>(QEDPartonIndices::SIGMAL);
						std::array s_dists{sigmaud, sigma, gluon, photon, sigmal};

						// for nf=4, we have an equal number of up/down quarks
						// so the numerator, Nup-Ndown = 0, and deltaNf = 0
						double deltaNf = 0;

						double cp = 0.5*((2./3.)*(2./3.) + (-1./3.)*(-1./3.));
						double cm = 0.5*((2./3.)*(2./3.) - (-1./3.)*(-1./3.));
						double fac_qcd = -2.0/_alpha_s.beta0();
						double fac_qed = -2.0/beta0qed;

						double L0QED = L1QED;
						double L0QCD = L1;

						for (uint s=1; s<_iterations; s++) {
							for (uint n=1; n<=s; n++) {
								auto pows = std::pow(L0QED,n)*std::pow(L0QCD,s-n)/factorial(n)/factorial(s-n);
								for (uint k=0; k<_grid.size()-1;k++) {
									double res1 = fac_qed*(
										cp*_grid.convolution(_A_QED(sigmaud,0,n-1), p0ff, k) +
										cm*_grid.convolution(_A_QED(sigma,0,n-1), p0ff, k) +
										0 +
										2*NC*_nf*(cp*deltaNf + cm)*_grid.convolution(_A_QED(photon,0,n-1), p0fy, k) +
										0
									);
									double res2 = fac_qed*(
										cm*_grid.convolution(_A_QED(sigmaud,0,n-1), p0ff, k) +
										cp*_grid.convolution(_A_QED(sigma,0,n-1), p0ff, k) +
										0 +
										2*NC*_nf*(cp + cm*deltaNf)*_grid.convolution(_A_QED(photon,0,n-1), p0fy, k) +
										0
									);
									double res4 = fac_qed*(
										cm*_grid.convolution(_A_QED(sigmaud,0,n-1), p0yf, k) +
										cp*_grid.convolution(_A_QED(sigma,0,n-1), p0yf, k) +
										0 +
										-3.0/4.0*beta0qed*_grid.convolution(_A_QED(photon,0,n-1), p0yy, k) +
										_grid.convolution(_A_QED(sigmal,0,n-1), p0yf, k)
									);
									double res5 = fac_qed*(
										0 +
										0 +
										0 +
										2.0*_nl*_grid.convolution(_A_QED(photon,0,n-1), p0fy, k) +
										_grid.convolution(_A_QED(sigmal,0,n-1), p0ff, k)
									);
									
									_S_QED(sigmaud,1,n,k) = res1;
									_S_QED(sigma,1,n,k) = res2;
									_S_QED(gluon,1,n,k) = 0; // obviously
									_S_QED(photon,1,n,k) = res4;
									_S_QED(sigmal,1,n,k) = res5;

									for (uint j : s_dists)
										resum_singlet[j][k] += _S_QED(j,1,n,k)*pows;
								}
							}

							uint n = 0;
							auto pows = std::pow(L0QED,n)*std::pow(L0QCD,s-n)/factorial(n)/factorial(s-n);
							for (uint k=0; k<_grid.size()-1;k++) {
								double res1 = fac_qcd*(
									_grid.convolution(_A_QED(sigmaud,0,n), p0ns, k) +
									deltaNf*(
										_grid.convolution(_A_QED(sigma,0,n), p0qq, k) +
										_grid.convolution(_A_QED(sigma,0,n), p0ns, k)) +
									deltaNf*_grid.convolution(_A_QED(gluon,0,n), p0qg, k) +
									0 +
									0
								);
								double res2 = fac_qcd*(
									0 +
									_grid.convolution(_A_QED(sigma,0,n), p0qq, k) +
									_grid.convolution(_A_QED(gluon,0,n), p0qg, k) +
									0 +
									0
								);
								double res3 = fac_qcd*(
									0 +
									_grid.convolution(_A_QED(sigma,0,n), p0gq, k) +
									_grid.convolution(_A_QED(gluon,0,n), p0gg, k) +
									0 +
									0
								);
								
								_S_QED(sigmaud,1,n,k) = res1;
								_S_QED(sigma,1,n,k) = res2;
								_S_QED(gluon,1,n,k) = res3;
								_S_QED(photon,1,n,k) = 0; // obviously
								_S_QED(sigmal,1,n,k) = 0; // obviously
								
								for (uint j : s_dists)
									resum_singlet[j][k] += _S_QED(j,1,n,k)*pows;
							}

							for (uint j : s_dists) {
								for (uint n=0; n<=s; ++n)
									std::ranges::copy(_S_QED(j,1,n), _S_QED(j,0,n).begin());
							}
						}
					}
					// non-singlet
					{
						std::array ns_dists{
							static_cast<uint>(QEDPartonIndices::UV),
							static_cast<uint>(QEDPartonIndices::DV),
							static_cast<uint>(QEDPartonIndices::SIGMADS),
							static_cast<uint>(QEDPartonIndices::SIGMAUC),
							static_cast<uint>(QEDPartonIndices::SIGMASB)
						};
						// aliases, TODO: change these
						double L0QED = L1QED;
						double L0QCD = L1;
						
						for (uint j : ns_dists) {
							auto& p0ns = getExpression(ExprName::P0ns);
							auto& p0ff     = getExpression(ExprName::P0ff);
							double fac_qcd = -2.0/_alpha_s.beta0();
							double fac_qed = -2.0/beta0qed;
		
							for (uint s=1; s<_iterations; s++) {
								for (uint n=1; n<=s; n++) {
									for (uint k=0; k<_grid.size()-1;k++) {
										_A_QED(j,1,n,k) = fac_qed*_grid.convolution(_A_QED(j,0,n-1), p0ff, k);
										resum_ns[j][k] += _A_QED(j,1,n,k)*std::pow(L0QED,n)*std::pow(L0QCD,s-n)/factorial(n)/factorial(s-n);
									}
								}
			
								for (uint k=0; k<_grid.size()-1;k++) {
									uint n = 0;
									_A_QED(j,1,n,k) = fac_qcd*_grid.convolution(_A_QED(j,0,n), p0ns, k);
									resum_ns[j][k] += _A_QED(j,1,n,k)
										*std::pow(L0QED,n)*std::pow(L0QCD,s-n)
										/factorial(n)/factorial(s-n);
								}
			
								for (uint n=0; n<=s; ++n)
									std::ranges::copy(_A_QED(j,1,n), _A_QED(j,0,n).begin());
							}
						}
					}
					fixDistributionsQED(resum_ns, resum_singlet, resum);
				} else {
					log(LOG_DEBUG, "DGLAP", "Starting singlet evolution and resummation...");
					evolveSinglet(resum_singlet, L1);
					log(LOG_DEBUG, "DGLAP", "Finished singlet evolution and resummation.");
					log(LOG_DEBUG, "DGLAP", "Starting non-singlet evolution and resummation...");
					_evol_type == EvolType::Exact ?
						evolveNonSinglet(resum_ns, L1, L2, L3, L4) :
						evolveNonSingletTrunc(resum_ns, L1);
					log(LOG_DEBUG, "DGLAP", "Finished non-singlet evolution and resummation.");

					log(LOG_DEBUG, "DGLAP", "Fixing distributions...");
					fixDistributions(resum_ns, resum_singlet, resum);
					log(LOG_DEBUG, "DGLAP", "Finished fixing distributions.");
				}

				// if we just resummed to a tabulated value,
				// _F contains our final distributions
				// we can just copy
				if (resum_tab) {
					log(LOG_DEBUG, "DGLAP", "Moving distributions into output array.");
					_F = std::move(resum);
				} else if (resum_threshold) {
					log(LOG_DEBUG, "DGLAP", "Moving distributions into the initial conditions of the next iteration.");
					// if we just resummed to a threshold energy,
					// then we need to recopy the resultant distributions
					// from the temporary array
					// back to the n=0 piece
					for (uint j=0; j<DISTS; ++j) {
						std::ranges::copy(resum[j], getNonSingletCoeffArray(j).begin());
					}
					for (uint j=0; j<=1; ++j)
						std::ranges::copy(resum[j*31], _S(0,j,0).begin());
				}
			} else { // if (alpha0 != alpha1)
				// if we've done no evolutions or anything,
				// we want to make sure we return correctly
				// the initial distributions
				if (last_loop && !performed_evolution) {
					fixDistributionsForce(resum);
					_F = std::move(resum);
					break;
				}
			}

			// +1 is the mass at the end of the current threshold,
			// so if +2 is zero then there is no next threshold,
			// meaning there is no use to do heavy flavor treatment/matching
			if (resum_threshold && _order>=2 && _alpha_s.masses(_nf+2)!=0.0)
				heavyFlavorTreatment();
			
		} // for (_nf=_nfi; ; _nf++)

		log(LOG_INFO, "DGLAP", "Done!");
		return _F;
	} // _evolve_function()


	std::vector<ArrayGrid> DGLAPSolver::calculateSubtractionPDFs()
	{
		if (_order < 2) {
			log(LOG_WARNING, "DGLAPSolver", "cannot create subtraction PDFs below NNLO");
			return {};
		}

		double as = _alpha1/(4.0*PI);
		double as2 = as*as;
		double as3 = as2*as;
		double mc = _initial_dist.masses(DIST_C);
		double mc2 = mc*mc;
		double mb = _initial_dist.masses(DIST_B);
		double mb2 = mb*mb;
		double qf = _Qf;
		double qf2 = _Qf*_Qf;
		

		auto zero_func = [](double,double){ return 0.0; };
	    auto a1qg_reg_func = [as](double lm, double nf, double x) {
			auto ome_reg = ome::AQg_reg[1];
			return ome_reg(lm, nf, x); };
		OpMatElemCustom a1hg(a1qg_reg_func, zero_func, zero_func);

		auto a2hq_reg_func = [as](double lm, double nf, double x) {
			auto ome_reg = ome::AQqPS_reg[2];
			return ome_reg(lm, nf, x); };
		OpMatElemCustom a2hq(a2hq_reg_func, zero_func, zero_func);
		
		auto a2hg_reg_func = [as](double lm, double nf, double x) {
			auto ome_reg = ome::AQg_reg[2];
			return ome_reg(lm, nf, x); };
		OpMatElemCustom a2hg(a2hg_reg_func, zero_func, zero_func);

		auto a3hq_reg_func = [as](double lm, double nf, double x) {
			auto ome_reg = ome::AQqPS_reg[3];
			return ome_reg(lm, nf, x); };
		OpMatElemCustom a3hq(a3hq_reg_func, zero_func, zero_func);
		
		auto a3hg_reg_func = [as](double lm, double nf, double x) {
			auto ome_reg = ome::AQg_reg[3];
			return ome_reg(lm, nf, x); };
		OpMatElemCustom a3hg(a3hg_reg_func, zero_func, zero_func);


		std::vector<ArrayGrid> subpdfs(20, ArrayGrid(_grid.size()));
		// charm
		if (mc2 >= qf2) {
			double L = std::log(mc2/qf2);
			OpMatElemN3LO::update(L, _nf);
			
			for (uint k=0; k<_grid.size(); ++k) {
				double ftilde1 = as*_grid.convolution(_F[0], a1hg, k);
				double ftilde2 = as2*(
					_grid.convolution(_F[31], a2hq, k) +
					_grid.convolution(_F[0], a2hg, k)
				);
				double ftilde3 = as3*(
					_grid.convolution(_F[31], a3hq, k) +
					_grid.convolution(_F[0], a3hg, k)
				);

				subpdfs[0][k] = ftilde1;
				subpdfs[1][k] = ftilde2;
				subpdfs[2][k] = ftilde3;
			
				subpdfs[3][k] = ftilde1 + ftilde2;
				subpdfs[4][k] = ftilde1 + ftilde2 + ftilde3;
			
				subpdfs[5][k] = _F[5][k] - ftilde1;
				subpdfs[6][k] = _F[5][k] - ftilde2;
				subpdfs[7][k] = _F[5][k] - ftilde3;
				
				subpdfs[8][k] = _F[5][k] - subpdfs[3][k];
				subpdfs[9][k] = _F[5][k] - subpdfs[4][k];
			}
		}
		// bottom
		if (mb2 >= qf2) {
		    double L = std::log(mb2/qf2);
			OpMatElemN3LO::update(L, _nf);
			
			for (uint k=0; k<_grid.size(); ++k) {
				double ftilde1 = as*_grid.convolution(_F[0], a1hg, k);
				double ftilde2 = as2*(
					_grid.convolution(_F[31], a2hq, k) +
					_grid.convolution(_F[0], a2hg, k)
				);
				double ftilde3 = as3*(
					_grid.convolution(_F[31], a3hq, k) +
					_grid.convolution(_F[0], a3hg, k)
				);

				subpdfs[10][k] = ftilde1;
				subpdfs[11][k] = ftilde2;
				subpdfs[12][k] = ftilde3;
			
				subpdfs[13][k] = ftilde1 + ftilde2;
				subpdfs[14][k] = ftilde1 + ftilde2 + ftilde3;
			
				subpdfs[15][k] = _F[5][k] - ftilde1;
				subpdfs[16][k] = _F[5][k] - ftilde2;
				subpdfs[17][k] = _F[5][k] - ftilde3;
				subpdfs[18][k]  = _F[5][k] - subpdfs[3][k];
				subpdfs[19][k]  = _F[5][k] - subpdfs[4][k];
			}
		}

		return subpdfs;
	}
	
} // namespace Candia2
