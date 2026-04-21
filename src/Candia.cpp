// Candia.cpp

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/ArrayGrid.hpp"
#include "Candia-v2/Grid.hpp"
#include "Candia-v2/SplittingFn.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"


// PDF indices
// 
// 0      gluons
// 1-6    quarks
// 7-12   antiquarks
// 13-18  q_i^-
// 19-24  q_i^+
// 25     q^(-)
// 26-30  q_{NS,1i}^(-)
// 31     q^(+)
// 32-36  q_{NS,1i}^(+)
// 37-39  

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
		uint order, Grid& grid, AlphaS const& alpha_s,
		double Qf, uint iterations, uint trunc_idx,
		Distribution const& initial_dist,
		double mur2_muf2) 
		: _order{order},  _grid{grid}, _Qf{Qf},
		  _alpha_s{alpha_s},
		  _mur2_muf2{mur2_muf2}, _log_mur2_muf2{std::log(mur2_muf2)}, _log_muf2_mur2{-_log_mur2_muf2},
		  _is_scale_difference{mur2_muf2 != 1.0},
		  _initial_dist{initial_dist},
		  _iterations{iterations}, _trunc_idx{trunc_idx}
	{
		log(::CANDIA_OPENING_TEXT);

		log(LOG_INFO, "DGLAP", "Evolving with log(mu_R / mu_F) = log({:.1}) = {:.4}.", _mur2_muf2, _log_mur2_muf2);

		switch(_order) {
			case 0: {
				if (_trunc_idx != 0) {
					_trunc_idx = 0; // LO has exact singlet solution, do not add additional terms
					log(LOG_WARNING, "DGLAP", "Specified value of the truncation index ({}) will be set to zero.", _trunc_idx);
				}

				_A = std::vector<std::vector<ArrayGrid>>{
					DISTS, std::vector<ArrayGrid>{
						2, ArrayGrid(_grid.size())
					}
				};
			} break;
			case 1: {
				_B = MultiDimArrayGrid_t<3>{
					DISTS, MultiDimArrayGrid_t<2>{
						2, MultiDimArrayGrid_t<1>{
							_iterations, ArrayGrid(_grid.size())
						}
					}
				};

			} break;
			case 2: {
				_C = MultiDimArrayGrid_t<4>{
					DISTS, MultiDimArrayGrid_t<3>{
						2, MultiDimArrayGrid_t<2>{
							_iterations, MultiDimArrayGrid_t<1>{
								_iterations, ArrayGrid(_grid.size())
							}
						}
					}
				};
			} break;
			case 3: {
				_D = MultiDimArrayGrid_t<5>{
					DISTS, MultiDimArrayGrid_t<4>{
						2, MultiDimArrayGrid_t<3>{
							_iterations, MultiDimArrayGrid_t<2>{
								_iterations, MultiDimArrayGrid_t<1>{
									_iterations, ArrayGrid(_grid.size())
								}
							}
						}
					}
				};

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
				log(LOG_DEBUG, "DGLAPSolver::DGLAPSolver()", "b  array: {}", func(_r1));
				log(LOG_DEBUG, "DGLAPSolver::DGLAPSolver()", "c  array: {}", func(_r1));
			} break;
			default: {
				log(LOG_INFO, "DGLAPSolver::DGLAPSolver()", "Found {} for the order, expected a value in range [0, 3].", _order);
			}
		}
		_S = decltype(_S){
			_trunc_idx+1, std::vector<std::vector<ArrayGrid>>{
				2, std::vector<ArrayGrid>{
					2, ArrayGrid(_grid.size())
				}
			}
		};
		log(LOG_INFO, "DGLAP", "Reserved space in coefficient arrays.");


		_F = std::vector<ArrayGrid>{
			DISTS, ArrayGrid(_grid.size())
		};

		setInitialConditions(initial_dist);
		log(LOG_INFO, "DGLAP", "Successfully filled coefficients with initial conditions.");
	}

	DGLAPSolver::~DGLAPSolver()
	{
		log(LOG_INFO, "DGLAP", "Exiting...");
		log(CANDIA_CLOSING_TEXT);
	}

	void DGLAPSolver::setInitialConditions(Distribution const& dist)
	{
		log(LOG_INFO, "DGLAP", "Setting initial conditions... ");

		dist.fillSingletCoeffs(
			[&](uint j, uint k) -> ArrayGrid::value_type& {
				return _S[0][j][0][k]; },
			_grid.points());
		dist.fillNonSingletCoeffs(
			[&](uint j, uint k) -> ArrayGrid::value_type& {
			switch (_order) {
				case 0: return _A[j][0][k]; break;
				case 1: return _B[j][0][0][k]; break;
				case 2: return _C[j][0][0][0][k]; break;
				case 3: return _D[j][0][0][0][0][k]; break;
				default:
					exit(EXIT_FAILURE); }},
			_grid.points());
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
			log(LOG_DEBUG, "DGLAP", "Loading Fortran versions of P2 splitting functions");
			createExpression<mvv_p2::P2nsm>(ExprName::P2nsm);
			createExpression<mvv_p2::P2nsp>(ExprName::P2nsp);
			createExpression<mvv_p2::P2nsv>(ExprName::P2nsv);
			createExpression<mvv_p2::P2qq>(ExprName::P2qq);
			createExpression<mvv_p2::P2qg>(ExprName::P2qg);
			createExpression<mvv_p2::P2gq>(ExprName::P2gq);
			createExpression<mvv_p2::P2gg>(ExprName::P2gg);
		} else {
			log(LOG_DEBUG, "DGLAP", "Loading C++ versions of P2 splitting functions");
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
		
		if (options.use_fortran_n3lo_splitfuncs) {
			log(LOG_DEBUG, "DGLAP", "Loading Fortran versions of P3 splitting functions");
			createExpression<mvv_p3::P3nsm>(ExprName::P3nsm);
			createExpression<mvv_p3::P3nsp>(ExprName::P3nsp);
			createExpression<mvv_p3::P3nsv>(ExprName::P3nsv);
			createExpression<mvv_p3::P3qq>(ExprName::P3qq);
			createExpression<mvv_p3::P3qg>(ExprName::P3qg);
			createExpression<mvv_p3::P3gq>(ExprName::P3gq);
			createExpression<mvv_p3::P3gg>(ExprName::P3gg);
		} else {
			log(LOG_DEBUG, "DGLAP", "Loading C++ versions of P3 splitting functions");
			createExpression<P3nsm>(ExprName::P3nsm);
			createExpression<P3nsp>(ExprName::P3nsp);
			createExpression<P3nsv>(ExprName::P3nsv);
			createExpression<P3qq>(ExprName::P3qq);
			createExpression<P3qg>(ExprName::P3qg);
			createExpression<P3gq>(ExprName::P3gq);
			createExpression<P3gg>(ExprName::P3gg);
		}

		uint imod = getOptions().n3lo_splitfunc_imod;
		SplittingFunction::setN3LOApproxType(imod);
		log(LOG_DEBUG, "DGLAP", "Setting N3LO approximation type (imod) = {}", imod);

		createExpression<OpMatElemN3LO>(ExprName::A3nsm, ome::AqqQNSEven);
		createExpression<OpMatElemN3LO>(ExprName::A3nsp, ome::AqqQNSOdd);
		createExpression<OpMatElemN3LO>(ExprName::A3gq, ome::AgqQ);
		createExpression<OpMatElemN3LO>(ExprName::A3gg, ome::AggQ);
		createExpression<OpMatElemN3LO>(ExprName::A3hq, ome::AQqPS);
		createExpression<OpMatElemN3LO>(ExprName::A3hg, ome::AQg);
		createExpression<OpMatElemN3LO>(ExprName::A3psqq, ome::AqqQPS);
		createExpression<OpMatElemN3LO>(ExprName::A3sqg, ome::AqgQ);
		createExpression<OpMatElemN3LO>(ExprName::A3PSshq, ome::AQqPSs);
    }

    void DGLAPSolver::setupCoefficients()
    {
		auto get_dist = [&](uint j, uint k) -> double& {
			if (options.use_truncated_nonsinglet_sol)
				return _S_NS[0][j][0][k];
			switch (_order) {
				case 0: return _A[j][0][k]; break;
				case 1: return _B[j][0][0][k]; break;
				case 2: return _C[j][0][0][0][k]; break;
				case 3: return _D[j][0][0][0][0][k]; break;
				default: throw "unreachable";
			}
		};

		for (uint k=0; k<_grid.size(); k++) {
			for (uint j=13; j<=18; j++)
				get_dist(j, k) = get_dist(j-12, k)-get_dist(j-6, k);
		
			get_dist(25, k)=0.;
			for (uint j=13; j<=18; j++)
				get_dist(25, k) += get_dist(j, k);

			for (uint j=26; j<=30; j++)
				get_dist(j, k) = get_dist(13, k)-get_dist(j-12, k);

			for (uint j=19; j<=24; j++)
				get_dist(j, k) = get_dist(j-18, k)+get_dist(j-12, k);

			_S[0][1][0][k] = 0.0;
			for (uint j=19; j<=24; j++)
				_S[0][1][0][k] += get_dist(j, k);
		
			for (uint j=32; j<=36; j++)
				get_dist(j, k)=get_dist(19, k)-get_dist(j-12, k);
		}
    }

	void DGLAPSolver::fixDistributions(
		std::vector<ArrayGrid>& resum_ns,
		std::vector<ArrayGrid>& resum_singlet,
		std::vector<ArrayGrid>& resum)
    {
        for (uint t=1; t<=_trunc_idx; ++t) {
			for (uint j=0; j<=1; ++j) {
				for (uint n=0; n<=1; ++n) 
					_S[t][j][n].zero();
			}
        }

		if (options.use_truncated_nonsinglet_sol) {
			for (uint t=1; t<=_trunc_idx; ++t) {
				for (uint j=0; j<_F.size(); ++j) {
					for (uint n=0; n<=1; ++n) 
						_S_NS[t][j][n].zero();
				}
			}
		}

		for (uint j=0; j<=1; ++j)
			resum[j*31] = resum_singlet[j*31];
		switch (_order) {
			case 0:
			case 1: {
				for (uint j=13; j<=12+_nf; j++)
					resum[j] = resum_ns[j];
				for (uint j=32; j<=30+_nf; j++)
					resum[j] = resum_ns[j];
			} break;
			case 2:
			case 3: {
				for (uint j=26; j<=24+_nf; ++j)
					resum[j] = resum_ns[j];
				for (uint j=32; j<=30+_nf; ++j)
					resum[j] = resum_ns[j];
				resum[25] = resum_ns[25];
			} break;
		}

        double Nf = static_cast<double>(_nf);
		for (uint k=0; k<_grid.size()-1;k++) {
			if (_order>=2) {
				resum[13][k]=resum[25][k];
				for (uint j=26; j<=24+_nf; j++)
					resum[13][k] += resum[j][k];
				resum[13][k] /= Nf;
				for (uint j=14; j<=12+_nf; j++)
					resum[j][k] = resum[13][k] - resum[j+12][k];
			}

			resum[19][k] = resum[31][k];
			for (uint j=32; j<=30+_nf; j++)
				resum[19][k] += resum[j][k];
			resum[19][k] /= Nf;

			for (uint j=20; j<=18+_nf; j++)
				resum[j][k] = resum[19][k] - resum[j+12][k];

			for (uint j=1; j<=_nf; j++) {
				resum[j][k]   =0.5*(resum[j+18][k] + resum[j+12][k]);
				resum[j+6][k] =0.5*(resum[j+18][k] - resum[j+12][k]);
			}

			if (_order<2) {
				resum[25][k]=0.0;
				for (uint j=13; j<=12+_nf; j++)
					resum[25][k] += resum[j][k];

				for (uint j=26; j<=24+_nf; j++)
					resum[j][k] = resum[13][k] - resum[j-12][k];
			}
		}
    }

	void DGLAPSolver::fixDistributionsForce(std::vector<ArrayGrid>& resum)
	{
		auto get_singlet_dist = [&](uint j) -> ArrayGrid& {
			if (!options.use_truncated_nonsinglet_sol)
				return _S[0][j][0];
			else
				return _S_NS[0][j][0];
		};
		auto get_nonsinglet_dist = [&](uint j) -> ArrayGrid& {
			if (options.use_truncated_nonsinglet_sol)
				return _S_NS[0][j][0];
			switch (_order) {
				case 0: return _A[j][0]; break;
				case 1: return _B[j][0][0]; break;
				case 2: return _C[j][0][0][0]; break;
				case 3: return _D[j][0][0][0][0]; break;
				default: throw "unreachable";
			}
		};


		for (uint j=0; j<=1; ++j)
			resum[j*31] = get_singlet_dist(j);
		switch (_order) {
			case 0:
			case 1: {
				for (uint j=13; j<=12+_nf; j++)
					resum[j] = get_nonsinglet_dist(j);
				for (uint j=32; j<=30+_nf; j++)
					resum[j] = get_nonsinglet_dist(j);
			} break;
			case 2:
			case 3: {
				for (uint j=26; j<=24+_nf; ++j)
					resum[j] = get_nonsinglet_dist(j);
				for (uint j=32; j<=30+_nf; ++j)
					resum[j] = get_nonsinglet_dist(j);
				resum[25] = get_nonsinglet_dist(25);
			} break;
		}

        double Nf = static_cast<double>(_nf);
		for (uint k=0; k<_grid.size()-1;k++) {
			if (_order>=2) {
				resum[13][k]=resum[25][k];
				for (uint j=26; j<=24+_nf; j++)
					resum[13][k] += resum[j][k];
				resum[13][k] /= Nf;
				for (uint j=14; j<=12+_nf; j++)
					resum[j][k] = resum[13][k] - resum[j+12][k];
			}

			resum[19][k] = resum[31][k];
			for (uint j=32; j<=30+_nf; j++)
				resum[19][k] += resum[j][k];
			resum[19][k] /= Nf;

			for (uint j=20; j<=18+_nf; j++)
				resum[j][k] = resum[19][k] - resum[j+12][k];

			for (uint j=1; j<=_nf; j++) {
				resum[j][k]   =0.5*(resum[j+18][k] + resum[j+12][k]);
				resum[j+6][k] =0.5*(resum[j+18][k] - resum[j+12][k]);
			}

			if (_order<2) {
				resum[25][k]=0.0;
				for (uint j=13; j<=12+_nf; j++)
					resum[25][k] += resum[j][k];

				for (uint j=26; j<=24+_nf; j++)
					resum[j][k] = resum[13][k] - resum[j-12][k];
			}
		}
	}

	void DGLAPSolver::setupTruncatedDistributions()
	{
		log(LOG_INFO, "Grid", "Using truncated ansatz for non-singlet sector.");
		
	    _S_NS = decltype(_S_NS){
			_trunc_idx+1, std::vector<std::vector<ArrayGrid>>{
				DISTS, std::vector<ArrayGrid>{
					2, ArrayGrid(_grid.size())
				}
			}
		};

		auto coeff_accessor = [&](uint j) -> ArrayGrid const& {
			switch (_order) {
				case 0: return _A[j][0];
				case 1: return _B[j][0][0];
				case 2: return _C[j][0][0][0];
				case 3: return _D[j][0][0][0][0];
				default: throw "unreachable";
			}
		};
		for (uint j=0; j<DISTS; ++j)
			_S_NS[0][j][0] = coeff_accessor(j);
		
		_A.clear();
		_B.clear();
		_C.clear();
		_D.clear();
	}

	std::vector<ArrayGrid> const& DGLAPSolver::evolve()
	{
	    log(LOG_INFO, "DGLAP", "Evolving to {} flavors.", _alpha_s.nff());
		using out_type = decltype(_F);
		loadAllExpressions();

		if (options.use_truncated_nonsinglet_sol)
			setupTruncatedDistributions();

		//std::array<double,1> Qtab{_Qf};
		out_type final_dists;

		// since we now only store two iterations at once,
		// we create these temporary arrays that will store the results of the resummation
		// since they were originally stored in the s=0 part that no longer exists
		// they are then copied to the _F final dists or the s=0 of the next iteration
		std::vector<ArrayGrid> resum_ns(DISTS, ArrayGrid(_grid.size()));
		std::vector<ArrayGrid> resum_singlet(DISTS, ArrayGrid(_grid.size()));
		std::vector<ArrayGrid> resum(DISTS, ArrayGrid(_grid.size()));
		
		bool performed_evolution = false;
		for (_nf=_alpha_s.nfi(); _nf<=_alpha_s.nff(); _nf++) {
			log(LOG_INFO, "DGLAP", "Setting nf={}", _nf);
			bool last_loop = _nf == _alpha_s.nff();

			log(LOG_INFO, "DGLAP", "Setting up distributions for evolution.");
			setupCoefficients();
			log(LOG_INFO, "DGLAP", "Finished setting up distributions for evolution.");

			// if the next mass is zero, we are already done
			if (_alpha_s.masses(_nf+1) == 0.0) {
				log(LOG_INFO, "DGLAP", "Next mass is zero. Quitting...");
				break;
			}

			// update all values
			_alpha_s.update(_nf);
			log(LOG_DEBUG, "DGLAP", "Loading relevant splitting function / OME values into cache");
			SplittingFunction::update(_nf, _alpha_s.beta0(), _log_muf2_mur2);
			OpMatElem::update(-_log_mur2_muf2, _nf);
			if (getOptions().cache_exprs) {
				log(LOG_INFO, "DGLAP", "Loading all expression values into caches...");
				_grid.useCachedExpressions();
				for (auto& expr : _expressions)
					expr->fill(_grid.points(), _grid.abscissae(), _grid.filler().getMappings(-1.0));
			}
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
			double r1 = _r1[_nf];
			double b = _b[_nf];
			double c = _c[_nf];
			double L1 = std::log(_alpha1/_alpha0);
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
					log(LOG_INFO, "DGLAP", "Encountered argument of square root that is <0. Doing analytic continuation.");
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

			log(LOG_DEBUG, "DGLAP::evolve()", "Values of log coeffs:");
			std::vector<std::pair<double, std::string_view>> coeffs{{L1, "L1"}, {L2, "L2"}, {L3, "L3"}, {L4, "L4"}};
			for (auto [x, xname] : coeffs)
				log(LOG_DEBUG, "DGLAP::evolve()", "  - {} = {: }", xname, x);
			
			log(LOG_INFO, "DGLAP", "Doing {} resummation", (resum_tab ? "tabulated" : "threshold" ));
			log(LOG_INFO, "DGLAP", "AlphaS: {} --> {}", _alpha0, _alpha1);

			// only do evolution if alphas are different
			// (i.e. energy scales are different)
			if (_alpha0 != _alpha1) {
				performed_evolution = true;
				log(LOG_INFO, "DGLAP", "Starting singlet evolution and resummation...");
#if ENABLE_THREADING
				evolveSingletThreaded(resum_singlet, L1);
#else
			    evolveSinglet(resum_singlet, L1);
#endif
				log(LOG_INFO, "DGLAP", "Finished singlet evolution and resummation.");

				log(LOG_INFO, "DGLAP", "Starting non-singlet evolution and resummation...");
#if ENABLE_THREADING
				options.use_truncated_nonsinglet_sol ?
					evolveNonSingletTruncThreaded(resum_ns, L1) :
					evolveNonSingletThreaded(resum_ns, L1, L2, L3, L4);
#else
				options.use_truncated_nonsinglet_sol ?
					evolveNonSingletTrunc(resum_ns, L1) :
					evolveNonSinglet(resum_ns, L1, L2, L3, L4);
#endif
				log(LOG_INFO, "DGLAP", "Finished non-singlet evolution and resummation.");

				log(LOG_INFO, "DGLAP", "Fixing distributions...");
				fixDistributions(resum_ns, resum_singlet, resum);
				log(LOG_INFO, "DGLAP", "Finished fixing distributions.");

				// if we just resummed to a tabulated value,
				// _F contains our final distributions
				// we can just copy
				if (resum_tab) {
					log(LOG_INFO, "DGLAP", "Moving distributions into output array.");
					_F = std::move(resum);
				} else if (resum_threshold) {
					log(LOG_INFO, "DGLAP", "Moving distributions into the initial conditions of the next iteration.");
					// if we just resummed to a threshold energy,
					// then we need to recopy the resultant distributions
					// from the temporary array
					// back to the n=0 piece
					for (uint j=0; j<DISTS; ++j) {
						if (options.use_truncated_nonsinglet_sol) {
							_S_NS[0][j][0] = resum[j];
						} else {
							switch (_order) {
								case 0: _A[j][0] 		   = resum[j]; break;
								case 1: _B[j][0][0] 	   = resum[j]; break;
								case 2: _C[j][0][0][0]     = resum[j]; break;
								case 3: _D[j][0][0][0][0]  = resum[j]; break;
							}
						}
					}
					for (uint j=0; j<=1; ++j)
						_S[0][j][0] = resum[j*31];
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
	} // evolve()


	std::vector<ArrayGrid> DGLAPSolver::calculateSubtractionPDFs()
	{
		if (_order < 2) {
			log(LOG_WARNING, "DGLAPSolver", "cannot create subtraction PDFs below NNLO");
			return {};
		}

		double as = _alpha1/4.0/PI;
		double as2 = as*as;
		double mb = _initial_dist.masses(DIST_B);
		double qf = _Qf;
		double L = std::log(std::pow(qf/mb, 2.0));
		
		OpMatElemN3LO::update(-L, _nf);

		auto plus_zero_func = [](double,double,double){ return 0.0; };
		auto delta_zero_func = [](double,double){ return 0.0; };
		
		// auto& p1qg = getExpression("P1qg");
	    auto a1qg_reg_func = [as](double lm, double nf, double x) {
			auto trunced = ome::AQg_reg.truncate(1);
			return trunced(as, lm, nf, x); };
		OpMatElemCustom a1hg(a1qg_reg_func, plus_zero_func, delta_zero_func);

		auto a2hq_reg_func = [as](double lm, double nf, double x) {
			auto trunced = ome::AQqPS_reg.truncate(2);
			return trunced(as, lm, nf, x); };
		OpMatElemCustom a2hq(a2hq_reg_func, plus_zero_func, delta_zero_func);
		
		auto a2hg_reg_func = [as](double lm, double nf, double x) {
			auto trunced = ome::AQg_reg.truncate(2);
			return trunced(as, lm, nf, x); };
		OpMatElemCustom a2hg(a2hg_reg_func, plus_zero_func, delta_zero_func);
		

		std::vector<ArrayGrid> subpdfs(4, ArrayGrid(_grid.size()));
		for (uint k=0; k<_grid.size(); ++k) {
			subpdfs[0][k] = as*_grid.convolution(_F[0], a1hg, k);
		    double ftilde2 = as2*(
				_grid.convolution(_F[31], a2hq, k) +
				_grid.convolution(_F[0], a2hg, k)
			);
			subpdfs[1][k] = subpdfs[0][k] + ftilde2;
			subpdfs[2][k] = std::abs(_F[5][k] - subpdfs[0][k]);
			subpdfs[3][k] = std::abs(_F[5][k] - subpdfs[1][k]);
		}

		return subpdfs;
	}

	void DGLAPSolverLHAPDF::evolve(
		double q0, double qf, double dq,
		DGLAPSolver::options_type const& dglap_options)
	{
	    double nsize = (qf-q0)/dq;
		int size = std::trunc(nsize) + 1;
		
		auto as_qs_view =
			std::views::iota(0, size)
			| std::views::transform([q0,dq](int x) -> double { return q0 + dq*x; });
		auto as_qs = std::vector<double>(as_qs_view.begin(), as_qs_view.end());
		as_qs.emplace_back(qf); // may be necessary
	    
		AlphaS alphas_all(_order, q0, qf, _dist.alpha0(), _mur2_muf2);
		alphas_all.setVFNS(_dist.masses(), _dist.nfi(), _dist.nff());
		std::vector<std::pair<double,double>> as_qvals = alphas_all.getValues(as_qs);
		std::vector<double> as_vals(as_qvals.size());
		std::ranges::transform(
			as_qvals, std::ranges::begin(as_vals),
			[](std::pair<double,double> const& p) -> double { return p.second; });
		_as_qs = std::move(as_qs);
		_as_vals = std::move(as_vals);

		_xvals = _grid.points();

		log(LOG_INFO, "DGLAPSolverLHAPDF", "Energy values: {}", vec_to_str(_as_qs));
		for (double q : _as_qs) {
			log(LOG_INFO, "DGLAPSolverLHAPDF", "Performing the evolution from {} to {}", q0, q);
			
			AlphaS alphas(_order, q0, q, _dist.alpha0(), _mur2_muf2);
			alphas.setVFNS(_dist.masses(), _dist.nfi(), _dist.nff());
			// alphas.setFFNS(4);

			DGLAPSolver solver(_order, _grid, alphas, q, _iterations, _trunc_idx, _dist, _mur2_muf2);
			solver.getOptions() = dglap_options;
			std::vector<ArrayGrid> F = solver.evolve();

			std::map<int, ArrayGrid> map{};
			map[-5] = F[11];
			map[-4] = F[10];
			map[-3] = F[9];
			map[-2] = F[8];
			map[-1] = F[7];
			map[1] = F[1];
			map[2] = F[2];
			map[3] = F[3];
			map[4] = F[4];
			map[5] = F[5];
			map[21] = F[0];

			std::vector<ArrayGrid> subtraction_pdfs = solver.calculateSubtractionPDFs();
			map[FTILDE1] = subtraction_pdfs[0];
			map[FTILDENLO] = subtraction_pdfs[1];
			map[DELTAF1] = subtraction_pdfs[2];
			map[DELTAFNLO] = subtraction_pdfs[3];

		    _all_pdfs.emplace_back(std::make_pair(q, map));
		}
		
		log(LOG_INFO, "DGLAPSolverLHAPDF", "Finished running.");
	}

	void DGLAPSolverLHAPDF::write()
	{
		log(LOG_INFO, "DGLAPSolverLHAPDF", "Spitting out the stuff...");
		std::filesystem::path pdfdir_path = std::filesystem::current_path()/_name;
		if (!std::filesystem::exists(pdfdir_path)) {
			if (!std::filesystem::create_directory(pdfdir_path))
				log(LOG_ERROR, "DGLAPSolverLHAPDF", "Failed to create the pdf outpout directory: '{}'", pdfdir_path.string());
		}
		
		std::string infofile_name = std::format("{}.info", _name);
		std::ifstream infofile_in(_infofile_in_path);
		std::string infofile_in_contents(
			std::istreambuf_iterator<char>(infofile_in),
			std::istreambuf_iterator<char>{});
		std::string pids_replace_str = "%PIDS%";
		std::string order_replace_str = "%ORDER%";
		std::string as_qs_replace_str = "%AS_QS%";
		std::string as_vals_replace_str = "%AS_VALS%";

		auto perform_replace = [](
			std::string& str,
			std::string const& what_to_replace,
			std::string const& replace_with)
		{
			std::string::size_type ipos;
			while ((ipos = str.find(what_to_replace)) != std::string::npos)
				str.replace(ipos, what_to_replace.size(), replace_with);
		};

		auto pids =
			_all_pdfs[0].second
			| std::views::transform([](std::pair<const int, ArrayGrid> const& p) -> int { return p.first; });

		perform_replace(infofile_in_contents, pids_replace_str, vec_to_str2(pids));
		perform_replace(infofile_in_contents, order_replace_str, std::to_string(_order));
		perform_replace(infofile_in_contents, as_qs_replace_str, vec_to_str2(_as_qs));
		perform_replace(infofile_in_contents, as_vals_replace_str, vec_to_str2(_as_vals));

		std::filesystem::path infofile_path = pdfdir_path/infofile_name;
		std::ofstream infofile(infofile_path);
		std::copy(infofile_in_contents.begin(), infofile_in_contents.end(), std::ostreambuf_iterator<char>(infofile));
		
		std::string datafile_name = std::format("{}_0000.dat", _name);
		std::filesystem::path datafile_path = pdfdir_path/datafile_name;
		std::ofstream datafile(datafile_path);
		datafile << "PdfType: central\n";
		datafile << "Format: lhagrid1\n";
		datafile << "---\n";
		datafile << "   " << vec_to_str2(_xvals, " ") << '\n';
		datafile << "   " << vec_to_str2(_as_qs, " ") << '\n';
		datafile << " " << vec_to_str2(pids, " ") << '\n';

		// the type is
	    // std::vector<std::pair<double, std::map<int,ArrayGrid>>>
		datafile << std::setprecision(10) << std::scientific;
		for (uint ix=0; ix<_xvals.size(); ++ix) {
			for (uint iq=0; iq<_as_qs.size(); ++iq) {
				datafile << "   ";
				for (int pid : pids)
					datafile << std::setw(17) << _all_pdfs[iq].second[pid][ix] << ' ';
				datafile << '\n';
			}
		}
		datafile << "---\n";
	}
	
} // namespace Candia2
