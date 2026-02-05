#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/ArrayGrid.hpp"
#include "Candia-v2/Grid.hpp"
#include "Candia-v2/SplittingFn.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"

#include <functional>
#include <memory>



// PDF indices
// 
// 0      gluons         g
// 1-6    quarks         u,d,s,c,b,t
// 7-12   antiquarks     au,ad,as,ac,ab,at
// 13-18  q_i^-          um,dm,sm,cm,bm,tm
// 19-24  q_i^+          up,dp,sp,cp,bp,tp
// 25     q^(-)
// 26-30  q_{NS,1i}^(-)  dd,sd,cd,bd,td
// 31     q^(+)
// 32-36  q_{NS,1i}^(+)  ds,ss,cs,bs,ts



namespace Candia2
{

	DGLAPSolver::DGLAPSolver(
		uint order, Grid const& grid, AlphaS const& alpha_s,
		double Qf, uint iterations, uint trunc_idx,
		Distribution const& initial_dist,
		double kr) 
		: _order{order},  _grid{grid}, _Qf{Qf},
		  _alpha_s{alpha_s},
		  _mur2_muf2{kr}, _log_mur2_muf2{std::log(kr)},
		  _iterations{iterations}, _trunc_idx{trunc_idx},
		  _use_n3lo_matching_conditions{true}
	{
		log(LOG_INFO, "DGLAP", "Evolving with log(mu_R / mu_F) = log({:.1}) = {:.4}.", _mur2_muf2, _log_mur2_muf2);

		switch(_order) {
			case 0: {
				_trunc_idx = 0; // LO has exact singlet solution, do not add additional terms

				_A = std::vector<std::vector<ArrayGrid>>{
					DISTS, std::vector<ArrayGrid>{
						2, ArrayGrid(grid.size())
					}
				};
			} break;
			case 1: {
				_B = MultiDimArrayGrid_t<3>{
					DISTS, MultiDimArrayGrid_t<2>{
						2, MultiDimArrayGrid_t<1>{
							_iterations, ArrayGrid(grid.size())
						}
					}
				};

			} break;
			case 2: {
				_C = MultiDimArrayGrid_t<4>{
					DISTS, MultiDimArrayGrid_t<3>{
						2, MultiDimArrayGrid_t<2>{
							_iterations, MultiDimArrayGrid_t<1>{
								_iterations, ArrayGrid(grid.size())
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
									_iterations, ArrayGrid(grid.size())
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
			} break;
			default: {
				log(LOG_INFO, "DGLAPSolver::DGLAPSolver()", "Found {} for the order, expected a value in range [0, 3].", order);
			}
		}
		_S = decltype(_S){
			trunc_idx+1, std::vector<std::vector<ArrayGrid>>{
				2, std::vector<ArrayGrid>{
					2, ArrayGrid(_grid.size())
				}
			}
		};
		log(LOG_INFO, "DGLAP", "Reserved space in coefficient arrays.");


		_F = std::vector<ArrayGrid>{
			DISTS, ArrayGrid(grid.size())
		};

		setInitialConditions(initial_dist);
	}

	DGLAPSolver::~DGLAPSolver()
	{
		log(LOG_INFO, "DGLAP", "Exiting...");
	}

	Expression& DGLAPSolver::getExpression(std::string_view name)
	{
		auto it = _expressions.find(name);
		if (it == _expressions.end())
			log(LOG_ERROR, "DGLAPSolver::getExpression()", "Expression '{}' does not exist.", name);
		return *it->second;
	}


	void DGLAPSolver::setInitialConditions(Distribution const& dist)
	{
		log(LOG_INFO, "DGLAP", "Setting initial conditions... ");

		/*
		for (uint k=0; k<_grid.size()-1; k++) {
			double x = _grid[k];
			_S[0][0][0][k] = dist.xg(x);
			_S[0][1][0][k] = dist.xqplus(x);
		}
		*/

		/*
		for (uint k=0; k<_grid.size()-1; k++) {
			double x = _grid[k];
			get_dist(1, k) = dist.xu(x);  // u
			get_dist(2, k) = dist.xd(x);  // d
			get_dist(3, k) = dist.xs(x);  // s
			get_dist(7, k) = dist.xub(x); // ub
			get_dist(8, k) = dist.xdb(x); // db
			get_dist(9, k) = dist.xs(x);  // sb ( = s)
		}
		*/

		dist.fillSingletCoeffs(
			[&](uint j, uint k) -> double& {
				return _S[0][j][0][k]; },
			_grid.points());
		dist.fillNonSingletCoeffs(
			[&](uint j, uint k) -> double& {
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
        createExpression<P0ns>("P0ns");
        createExpression<P0qq>("P0qq");
        createExpression<P0qg>("P0qg");
        createExpression<P0gq>("P0gq");
        createExpression<P0gg>("P0gg");
    
        if (_order >= 1) {
		    createExpression<P1nsm>("P1nsm");
            createExpression<P1nsp>("P1nsp");
            createExpression<P1qq>("P1qq");
            createExpression<P1qg>("P1qg");
            createExpression<P1gq>("P1gq");
            createExpression<P1gg>("P1gg");
        }
        if (_order >= 2) {
            createExpression<P2nsm>("P2nsm");
            createExpression<P2nsp>("P2nsp");
            createExpression<P2nsv>("P2nsv");
            createExpression<P2qq>("P2qq");
            createExpression<P2qg>("P2qg");
            createExpression<P2gq>("P2gq");
            createExpression<P2gg>("P2gg");
			
            createExpression<A2ns>("A2ns");
            createExpression<A2gq>("A2gq");
            createExpression<A2gg>("A2gg");
            createExpression<A2hq>("A2hq");
            createExpression<A2hg>("A2hg");
        }
        if (_order >= 3) {
            createExpression<P3nsm>("P3nsm");
            createExpression<P3nsp>("P3nsp");
            createExpression<P3nsv>("P3nsv");
            createExpression<P3qq>("P3qq");
            createExpression<P3qg>("P3qg");
            createExpression<P3gq>("P3gq");
            createExpression<P3gg>("P3gg");

            createExpression<OpMatElemN3LO>("A3nsm", ome::AqqQNSEven);
            createExpression<OpMatElemN3LO>("A3nsp", ome::AqqQNSOdd);
            createExpression<OpMatElemN3LO>("A3gq", ome::AgqQ);
            createExpression<OpMatElemN3LO>("A3gg", ome::AggQ);
            createExpression<OpMatElemN3LO>("A3hq", ome::AQqPS);
            createExpression<OpMatElemN3LO>("A3hg", ome::AQg);
            createExpression<OpMatElemN3LO>("A3psqq", ome::AqqQPS);
            createExpression<OpMatElemN3LO>("A3sqg", ome::AqgQ);
        }
        
    }

    void DGLAPSolver::setupCoefficients()
    {
		auto get_dist = [&](uint j, uint k) -> double& {
			switch (_order) {
				case 0: return _A[j][0][k]; break;
				case 1: return _B[j][0][0][k]; break;
				case 2: return _C[j][0][0][0][k]; break;
				case 3: return _D[j][0][0][0][0][k]; break;
				default:
					exit(EXIT_FAILURE); // unreachable
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

	void DGLAPSolver::fixDistributions(bool resum_tab, bool resum_threshold, std::vector<ArrayGrid>& temp_arr, std::vector<ArrayGrid>& temp_arr_singlet)
    {
        for (uint t=1; t<=_trunc_idx; ++t) {
            for (uint j=0; j<=1; ++j) {
                for (uint n=0; n<=1; ++n)
                    _S[t][j][n].zero();
            }
        }

        double Nf = static_cast<double>(_nf);
        if (resum_tab) {
            for (uint k=0; k<_grid.size()-1;k++) {
                if (_order>=2) {
                    _F[13][k]=_F[25][k];
                    for (uint j=26; j<=24+_nf; j++)
                        _F[13][k] += _F[j][k];
                    _F[13][k] /= Nf;
                    for (uint j=14; j<=12+_nf; j++)
                        _F[j][k] = _F[13][k] - _F[j+12][k];
                }

                _F[19][k] = _F[31][k];
                for (uint j=32; j<=30+_nf; j++)
                    _F[19][k] += _F[j][k];
                _F[19][k] /= Nf;

                for (uint j=20; j<=18+_nf; j++)
                    _F[j][k] = _F[19][k] - _F[j+12][k];

                for (uint j=1; j<=_nf; j++) {
                    _F[j][k]   =0.5*(_F[j+18][k] + _F[j+12][k]);
                    _F[j+6][k] =0.5*(_F[j+18][k] - _F[j+12][k]);
                }

                if (_order<2) {
                    _F[25][k]=0.0;
                    for (uint j=13; j<=12+_nf; j++)
                        _F[25][k] += _F[j][k];

                    for (uint j=26; j<=24+_nf; j++)
                        _F[j][k] = _F[13][k] - _F[j-12][k];
                }
            }
        }
        else if (resum_threshold)
        {
            switch (_order) {
                case 0: {
					for (uint k=0; k<_grid.size()-1; k++) {
                        temp_arr[19][k] = temp_arr_singlet[31][k];
                        for (uint j=32; j<=30+_nf; j++)
                            temp_arr[19][k] += temp_arr[j][k];
                        temp_arr[19][k] /= Nf;

                        for (uint j=20; j<=18+_nf; j++)
                            temp_arr[j][k] = temp_arr[19][k] - temp_arr[j+12][k];

                        for (uint j=1; j<=_nf; j++) {
                            temp_arr[j][k]   = 0.5*(temp_arr[j+18][k] + temp_arr[j+12][k]);
                            temp_arr[j+6][k] = 0.5*(temp_arr[j+18][k] - temp_arr[j+12][k]);
                        }
                    }
				}; break;
                case 1: {
                    for (uint k=0; k<_grid.size()-1; k++) {
                        temp_arr[19][k] = temp_arr_singlet[31][k];
                        for (uint j=32; j<=30+_nf; j++)
                            temp_arr[19][k] += temp_arr[j][k];
                        temp_arr[19][k] /= Nf;

                        for (uint j=20; j<=18+_nf; j++)
                            temp_arr[j][k] = temp_arr[19][k] - temp_arr[j+12][k];

                        for (uint j=1; j<=_nf; j++) {
                            temp_arr[j][k]   = 0.5*(temp_arr[j+18][k] + temp_arr[j+12][k]);
                            temp_arr[j+6][k] = 0.5*(temp_arr[j+18][k] - temp_arr[j+12][k]);
                        }
                    }
                }; break;
                case 2: {
					for (uint k=0; k<_grid.size()-1; k++) {
                        temp_arr[13][k] = temp_arr[25][k];
                        for (uint j=26; j<=24+_nf; j++)
                            temp_arr[13][k] += temp_arr[j][k];
                        temp_arr[13][k] /= Nf;

                        for (uint j=14;j<=12+_nf;j++)
                            temp_arr[j][k] = temp_arr[13][k] - temp_arr[j+12][k];

                        temp_arr[19][k] = temp_arr_singlet[31][k];
                        for (uint j=32; j<=30+_nf; j++)
                            temp_arr[19][k] += temp_arr[j][k];
                        temp_arr[19][k] /= Nf;

                        for (uint j=20; j<=18+_nf; j++)
                            temp_arr[j][k] = temp_arr[19][k] - temp_arr[j+12][k];

                        for (uint j=1; j<=_nf; j++) {
                            temp_arr[j][k]  =0.5*(temp_arr[j+18][k]+temp_arr[j+12][k]);
                            temp_arr[j+6][k]=0.5*(temp_arr[j+18][k]-temp_arr[j+12][k]);
                        }
                    }
				}; break;
                case 3: {
                    for (uint k=0; k<_grid.size()-1; k++) {
                        temp_arr[13][k] = temp_arr[25][k];
                        for (uint j=26; j<=24+_nf; j++)
                            temp_arr[13][k] += temp_arr[j][k];
                        temp_arr[13][k] /= Nf;

                        for (uint j=14;j<=12+_nf;j++)
                            temp_arr[j][k] = temp_arr[13][k] - temp_arr[j+12][k];

                        temp_arr[19][k] = temp_arr_singlet[31][k];
                        for (uint j=32; j<=30+_nf; j++)
                            temp_arr[19][k] += temp_arr[j][k];
                        temp_arr[19][k] /= Nf;

                        for (uint j=20; j<=18+_nf; j++)
                            temp_arr[j][k] = temp_arr[19][k] - temp_arr[j+12][k];

                        for (uint j=1; j<=_nf; j++) {
                            temp_arr[j][k]  =0.5*(temp_arr[j+18][k]+temp_arr[j+12][k]);
                            temp_arr[j+6][k]=0.5*(temp_arr[j+18][k]-temp_arr[j+12][k]);
                        }
                    }
                } break;
            }
        }
    }

	auto DGLAPSolver::evolve() -> decltype(_F)
	{
		log(LOG_INFO, "DGLAP", "Evolving to {} flavors.", _alpha_s.nff());
		using out_type = decltype(_F);
		loadAllExpressions();

		//std::array<double,1> Qtab{_Qf};
		out_type final_dists;

		// temp array for the threshold summation
		// originally, we wrote directly to the n=0 piece
		// of the coefficient array,
		// but now we only have the two pieces that continually
		// evolve forward
		// so we stick everything into this temp array
		// during the evolution, then move it into the n=0
		// part after the full evolution
		std::vector<ArrayGrid> temp_arr(DISTS, ArrayGrid(_grid.size()));
		std::vector<ArrayGrid> temp_arr_singlet(DISTS, ArrayGrid(_grid.size()));
			
		// since the only difference during the evolution/resummation to
		// the tabulated energy or the threshold is what array we append to, 
		// I create two reference arrays (one for singlet one for non-singlet)
		// that are updated depending on what we are evolving to
		// that way we don't have to evaluate if's inside the grid loop
		std::reference_wrapper<std::vector<ArrayGrid>> arr{std::ref(temp_arr)};
		std::reference_wrapper<std::vector<ArrayGrid>> arr_singlet{std::ref(temp_arr_singlet)};

		for (_nf=_alpha_s.nfi(); _nf<=_alpha_s.nff(); _nf++) {
			log(LOG_INFO, "DGLAP", "Setting nf={}", _nf);

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
			SplittingFunction::update(_nf, _alpha_s.beta0());
			for (auto& [_, expr] : _expressions) {
				if (_grid.splitIntervals())
					expr->fill(_grid.points(), _grid.allAbscissae());
				else
					expr->fill(_grid.points(), _grid.abscissae());
			}
			_alpha0 = _alpha_s.post(_nf);
			_alpha1 = _alpha_s.pre(_nf+1);
			bool resum_tab = _alpha_s.resumTabulated();
			bool resum_threshold = !resum_tab;
			
			// alpha1 needs to be manually calculated
			// if we are evolving to a tabulated energy
			if (resum_tab) {
				_alpha1 = _alpha_s.evaluate(_alpha_s.masses(_nf), _Qf, _alpha0);
				arr = std::ref(_F);
				arr_singlet = std::ref(_F);
			} else {
				arr = std::ref(temp_arr);
				arr_singlet = std::ref(temp_arr_singlet);
			}

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
			
			log(LOG_INFO, "DGLAP", "Doing {} resummation", (resum_tab ? "tabulated" : "threshold" ));
			log(LOG_INFO, "DGLAP", "AlphaS: {} --> {}", _alpha0, _alpha1);

			// only do evolution if alphas are different
			// (i.e. energy scales are different)
			if (_alpha0 != _alpha1) {
				log(LOG_INFO, "DGLAP", "Starting singlet evolution and resummation...");
#if ENABLE_THREADING
				evolveSingletThreaded(arr_singlet, L1);
#else
				evolveSinglet(arr_singlet, L1);
#endif
				log(LOG_INFO, "DGLAP", "Finished singlet evolution and resummation.");

				log(LOG_INFO, "DGLAP", "Starting non-singlet evolution and resummation...");
#if ENABLE_THREADING
				evolveNonSingletThreaded(arr, L1, L2, L3, L4);
#else
				evolveNonSinglet(arr, L1, L2, L3, L4);
#endif
				log(LOG_INFO, "DGLAP", "Finished non-singlet evolution and resummation.");

				log(LOG_INFO, "DGLAP", "Fixing distributions...");
				fixDistributions(resum_tab, resum_threshold, temp_arr, temp_arr_singlet);
				log(LOG_INFO, "DGLAP", "Finished fixing distributions.");

				// if we just resummed to a tabulated value,
				// _F contains our final distributions
				// we can just copy
				if (resum_tab) {
					final_dists = _F;
				} else if (resum_threshold) {
					// if we just resummed to a threshold energy,
					// then we need to recopy the resultant distributions
					// from the temporary array
					// back to the n=0 piece
					for (uint j=0; j<DISTS; ++j) {
						switch (_order) {
							case 0: _A[j][0] 		   = temp_arr[j]; break;
							case 1: _B[j][0][0] 	   = temp_arr[j]; break;
							case 2: _C[j][0][0][0]    = temp_arr[j]; break;
							case 3: _D[j][0][0][0][0] = temp_arr[j]; break;
						}
						
					}
					for (uint j=0; j<=1; ++j)
						_S[0][j][0] = temp_arr_singlet[j*31];
				}
			} // if (alpha0 != alpha1)

			if (resum_threshold && _order>=2 && _alpha_s.masses(_nf+2)!=0.0)
				heavyFlavorTreatment();
			
		} // for (_nf=_nfi; ; _nf++)

		log(LOG_INFO, "DGLAP", "Done!");
		return final_dists;
	} // evolve()

	
} // namespace Candia2
