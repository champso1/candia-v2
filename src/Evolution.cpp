// EvolutionThread.cpp

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/ArrayGrid.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/Math.hpp"

#include <cmath>
#include <functional>
#include <thread>
#include <execution>

namespace Candia2
{
	thread_local int thread_index = -1;
	
	void DGLAPSolver::evolveSinglet(std::reference_wrapper<std::vector<ArrayGrid>> arr, double L1)
	{
        for (uint j=0; j<=1; ++j)
			std::ranges::copy(_S(0,j,0), arr.get()[j*31].begin());

		auto grid_idxs = std::ranges::views::iota(uint{0}, _grid.size()-1);

		startLogIterations();
        switch (_order) {
            case 0: {
				auto& p0qq = getExpression(ExprName::P0qq);
				auto& p0qg = getExpression(ExprName::P0qg);
				auto& p0gq = getExpression(ExprName::P0gq);
				auto& p0gg = getExpression(ExprName::P0gg);
				
                for (uint n=1; n<_iterations; n++) {
					logIterations(n, _iterations-1, "SingletLO");

					double fac = std::pow(L1, n)/factorial(n);

					std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k) {
						_S(0,1,1,k) =
							recrelS_1(_S(0,1,0), k, p0qq) +
							recrelS_1(_S(0,0,0), k, p0qg);
                        _S(0,0,1,k) =
							recrelS_1(_S(0,1,0), k, p0gq) +
							recrelS_1(_S(0,0,0), k, p0gg);
					});

					for (uint j=0; j<=1; j++) {
						for (uint k=0; k<_grid.size()-1; k++)
							arr.get()[j*31][k] += _S(0,j,1,k) * fac;
						std::ranges::copy(_S(0,j,1), _S(0,j,0).begin());
					}
                }
            } break;
            case 1: {
				auto& p0qq = getExpression(ExprName::P0qq);
				auto& p0qg = getExpression(ExprName::P0qg);
				auto& p0gq = getExpression(ExprName::P0gq);
				auto& p0gg = getExpression(ExprName::P0gg);
				auto& p1qq = getExpression(ExprName::P1qq);
				auto& p1qg = getExpression(ExprName::P1qg);
				auto& p1gq = getExpression(ExprName::P1gq);
				auto& p1gg = getExpression(ExprName::P1gg);
				
			    for (uint n=1; n<_iterations; n++) {
					// log(LOG_INFO, "DGLAP", "NLO Singlet Iteration {}", n);
					logIterations(n, _iterations-1, "SingletNLO");

					double fac = std::pow(L1, n)/factorial(n);
					
                    // LO piece (non truncated)
                    std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
                        _S(0,1,1,k) = 
							recrelS_1(_S(0,1,0), k, p0qq) +
							recrelS_1(_S(0,0,0), k, p0qg);
                        _S(0,0,1,k) = 
							recrelS_1(_S(0,1,0), k, p0gq) +
							recrelS_1(_S(0,0,0), k, p0gg);
                    });

                    // new NLO piece non-convolution
                    for (uint j=0; j<=1; j++) {
                        std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
							_S(1,j,1,k) = -_S(0,j,1,k) * _alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()) - _S(1,j,0,k);
						});     
                    }

					std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
						auto fac1=recrelS_2(_S(1,1,0), _S(0,1,0), k, p0qq, p1qq);
						auto fac2=recrelS_2(_S(1,0,0), _S(0,0,0), k, p0qg, p1qg);
						_S(1,1,1,k) +=
						    fac1 +
						    fac2;
						_S(1,0,1,k) +=
							recrelS_2(_S(1,1,0), _S(0,1,0), k, p0gq, p1gq) +
							recrelS_2(_S(1,0,0), _S(0,0,0), k, p0gg, p1gg);
					});

                    // NLO truncation terms
                    for (uint t=2; t<=_trunc_idx; ++t) {
                        double T = static_cast<double>(t);

                        // non-convolution piece:
                        for (uint j=0; j<=1; j++) {
							std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
								_S(t,j,1,k) =
                                    - (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S(t-1,j,1,k)
                                    - T*_S(t,j,0,k)
                                    - (T-1.0)*(_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S(t-1,j,0,k);
							});
                        }

						std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
							_S(t,1,1,k) += 
								recrelS_2(_S(t,1,0), _S(t-1,1,0), k, p0qq, p1qq) +
								recrelS_2(_S(t,0,0), _S(t-1,0,0), k, p0qg, p1qg);
							_S(t,0,1,k) += 
								recrelS_2(_S(t,1,0), _S(t-1,1,0), k, p0gq, p1gq) +
								recrelS_2(_S(t,0,0), _S(t-1,0,0), k, p0gg, p1gg);
						});
                    }

                    for (uint t=0; t<=_trunc_idx; ++t) {
						double as_fac = std::pow(_alpha1, t);
                        for (uint j=0; j<=1; ++j) {
                            std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
                                arr.get()[j*31][k] += _S(t,j,1,k) * as_fac * fac;
							});
							std::ranges::copy(_S(t,j,1), _S(t,j,0).begin());
                        }
                    }
                }
            } break;
            case 2: {
				auto& p0qq = getExpression(ExprName::P0qq);
				auto& p0qg = getExpression(ExprName::P0qg);
				auto& p0gq = getExpression(ExprName::P0gq);
				auto& p0gg = getExpression(ExprName::P0gg);
				auto& p1qq = getExpression(ExprName::P1qq);
				auto& p1qg = getExpression(ExprName::P1qg);
				auto& p1gq = getExpression(ExprName::P1gq);
				auto& p1gg = getExpression(ExprName::P1gg);
				auto& p2qq = getExpression(ExprName::P2qq);
				auto& p2qg = getExpression(ExprName::P2qg);
				auto& p2gq = getExpression(ExprName::P2gq);
				auto& p2gg = getExpression(ExprName::P2gg);
				
                for (uint n=1; n<_iterations; n++) {
					// log(LOG_INFO, "DGLAP", "NNLO Singlet Iteration {}", n);
					logIterations(n, _iterations-1, "SingletNNLO");

                    // LO piece (non truncated)
                    for (uint k=0; k<_grid.size()-1; k++) {
                        _S(0,1,1,k) =
							recrelS_1(_S(0,1,0), k, p0qq) +
							recrelS_1(_S(0,0,0), k, p0qg);
                        _S(0,0,1,k) =
							recrelS_1(_S(0,1,0), k, p0gq) +
							recrelS_1(_S(0,0,0), k, p0gg);
                    }

                    // new NLO piece non-convolution
                    for (uint j=0; j<=1; j++) {
                        for (uint k=0; k<_grid.size()-1; k++) 
                            _S(1,j,1,k) = -_S(0,j,1,k) * _alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()) - _S(1,j,0,k);
                    }

					std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
						_S(1,1,1,k) += 
							recrelS_2(_S(1,1,0), _S(0,1,0), k, p0qq, p1qq) +
							recrelS_2(_S(1,0,0), _S(0,0,0), k, p0qg, p1qg);
						_S(1,0,1,k) += 
							recrelS_2(_S(1,1,0), _S(0,1,0), k, p0gq, p1gq) +
							recrelS_2(_S(1,0,0), _S(0,0,0), k, p0gg, p1gg);
					});

                    // new NNLO piece non-convolution
                    for (uint j=0; j<=1; j++) {
                        for (uint k=0; k<_grid.size()-1; k++) {
                            _S(2,j,1,k) =
                                - (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S(1,j,1,k)
                                - (_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S(0,j,1,k)
                                - 2.0*_S(2,j,0,k)
                                - (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S(1,j,0,k);
						}
                    }

					std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
						_S(2,1,1,k) += 
							recrelS_3(_S(2,1,0), _S(1,1,0), _S(0,1,0), k, p0qq, p1qq, p2qq) +
							recrelS_3(_S(2,0,0), _S(1,0,0), _S(0,0,0), k, p0qg, p1qg, p2qg);
						_S(2,0,1,k) += 
							recrelS_3(_S(2,1,0), _S(1,1,0), _S(0,1,0), k, p0gq, p1gq, p2gq) +
							recrelS_3(_S(2,0,0), _S(1,0,0), _S(0,0,0), k, p0gg, p1gg, p2gg);
					});

                    for (uint t=3; t<=_trunc_idx; ++t) {
                        double T = static_cast<double>(t);

                        // non-convolution piece:
                        for (uint j=0; j<=1; j++) {
                            for (uint k=0; k<_grid.size()-1; k++) {
                                _S(t,j,1,k) =
                                    - (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S(t-1,j,1,k)
                                    - (_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S(t-2,j,1,k)
                                    - T*_S(t,j,0,k)
                                    - (T-1.0)*(_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S(t-1,j,0,k)
                                    - (T-2.0)*(_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S(t-2,j,0,k);
							}
                        }

						std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
							_S(t,1,1,k) += 
								recrelS_3(_S(t,1,0), _S(t-1,1,0), _S(t-2,1,0), k, p0qq, p1qq, p2qq) +
								recrelS_3(_S(t,0,0), _S(t-1,0,0), _S(t-2,0,0), k, p0qg, p1qg, p2qg);
							_S(t,0,1,k) += 
								recrelS_3(_S(t,1,0), _S(t-1,1,0), _S(t-2,1,0), k, p0gq, p1gq, p2gq) +
								recrelS_3(_S(t,0,0), _S(t-1,0,0), _S(t-2,0,0), k, p0gg, p1gg, p2gg);
						});
                    }

                    for (uint t=0; t<=_trunc_idx; ++t) {
                        for (uint j=0; j<=1; ++j) {
                            for (uint k=0; k<_grid.size()-1; k++)
                                arr.get()[j*31][k] += _S(t,j,1,k) * std::pow(_alpha1, t) * std::pow(L1, n)/factorial(n);
							std::ranges::copy(_S(t,j,1), _S(t,j,0).begin());
                        }
                    }
                }
            } break;
            case 3: {
				auto& p0qq = getExpression(ExprName::P0qq);
				auto& p0qg = getExpression(ExprName::P0qg);
				auto& p0gq = getExpression(ExprName::P0gq);
				auto& p0gg = getExpression(ExprName::P0gg);
				auto& p1qq = getExpression(ExprName::P1qq);
				auto& p1qg = getExpression(ExprName::P1qg);
				auto& p1gq = getExpression(ExprName::P1gq);
				auto& p1gg = getExpression(ExprName::P1gg);
				auto& p2qq = getExpression(ExprName::P2qq);
				auto& p2qg = getExpression(ExprName::P2qg);
				auto& p2gq = getExpression(ExprName::P2gq);
				auto& p2gg = getExpression(ExprName::P2gg);
				auto& p3qq = getExpression(ExprName::P3qq);
				auto& p3qg = getExpression(ExprName::P3qg);
				auto& p3gq = getExpression(ExprName::P3gq);
				auto& p3gg = getExpression(ExprName::P3gg);

				auto beta0 = _alpha_s.beta0();
				auto beta1 = _alpha_s.beta1();
				auto beta2 = _alpha_s.beta2();
				auto beta3 = _alpha_s.beta3();
				
				auto fac1 = beta1/(4.0*PI*beta0);
				auto fac2 = beta2/(16.0*PI_2*beta0);
				auto fac3 = beta3/(64.0*PI_3*beta0);
				
			    for (uint n=1; n<_iterations; n++) {
					// log(LOG_INFO, "DGLAP", "N3LO Singlet Iteration {}", n);
					logIterations(n, _iterations-1, "SingletN3LO");

                    // LO piece
                    for (uint k=0; k<_grid.size()-1; k++) {
                        _S(0,1,1,k) =
							recrelS_1(_S(0,1,0), k, p0qq) +
							recrelS_1(_S(0,0,0), k, p0qg);
                        _S(0,0,1,k) =
							recrelS_1(_S(0,1,0), k, p0gq) +
							recrelS_1(_S(0,0,0), k, p0gg);
                    }

                    // new NLO piece non-convolution
                    for (uint j=0; j<=1; j++) {
                        for (uint k=0; k<_grid.size()-1; k++)
                            _S(1,j,1,k) = -_S(0,j,1,k)*fac1 - _S(1,j,0,k);
                    }

                    // new NLO piece convolution
                    std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
						_S(1,1,1,k) += 
							recrelS_2(_S(1,1,0), _S(0,1,0), k, p0qq, p1qq) +
							recrelS_2(_S(1,0,0), _S(0,0,0), k, p0qg, p1qg);
						_S(1,0,1,k) += 
							recrelS_2(_S(1,1,0), _S(0,1,0), k, p0gq, p1gq) +
							recrelS_2(_S(1,0,0), _S(0,0,0), k, p0gg, p1gg);
					});

                    // new NNLO piece non-convolution
                    for (uint j=0; j<=1; j++) {
                        for (uint k=0; k<_grid.size()-1; k++) {
                            _S(2,j,1,k) =
                                - fac1*_S(1,j,1,k)
                                - fac2*_S(0,j,1,k)
                                - 2.0*_S(2,j,0,k)
                                - fac1*_S(1,j,0,k);
						}
                    }

                    // new NNLO piece non-convolution
                    std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
						_S(2,1,1,k) += 
							recrelS_3(_S(2,1,0), _S(1,1,0), _S(0,1,0), k, p0qq, p1qq, p2qq) +
							recrelS_3(_S(2,0,0), _S(1,0,0), _S(0,0,0), k, p0qg, p1qg, p2qg);
						_S(2,0,1,k) += 
							recrelS_3(_S(2,1,0), _S(1,1,0), _S(0,1,0), k, p0gq, p1gq, p2gq) +
							recrelS_3(_S(2,0,0), _S(1,0,0), _S(0,0,0), k, p0gg, p1gg, p2gg);
					});

                    // new N3LO piece non-convolution
                    for (uint j=0; j<=1; j++) {
                        for (uint k=0; k<_grid.size()-1; k++) {
                            _S(3,j,1,k) =
                                - fac1*_S(2,j,1,k)
                                - fac2*_S(1,j,1,k)
                                - fac3*_S(0,j,1,k)
                                - 3.0*_S(3,j,0,k)
                                - 2.0*fac1*_S(2,j,0,k)
                                - fac2*_S(1,j,0,k);
						}
                    }
					
					std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
						_S(3,1,1,k) += 
							recrelS_4(_S(3,1,0), _S(2,1,0), _S(1,1,0), _S(0,1,0), k, p0qq, p1qq, p2qq, p3qq) +
							recrelS_4(_S(3,0,0), _S(2,0,0), _S(1,0,0), _S(0,0,0), k, p0qg, p1qg, p2qg, p3qg);
						_S(3,0,1,k) += 
							recrelS_4(_S(3,1,0), _S(2,1,0), _S(1,1,0), _S(0,1,0), k, p0gq, p1gq, p2gq, p3gq) +
							recrelS_4(_S(3,0,0), _S(2,0,0), _S(1,0,0), _S(0,0,0), k, p0gg, p1gg, p2gg, p3gg);
					});

                    for (uint t=4; t<=_trunc_idx; ++t) {
                        double T = static_cast<double>(t);

                        // non-convolution piece:
                        for (uint j=0; j<=1; j++) {
                            for (uint k=0; k<_grid.size()-1; k++) {
                                _S(t,j,1,k) =
                                    - fac1*_S(t-1,j,1,k)
                                    - fac2*_S(t-2,j,1,k)
                                    - fac3*_S(t-3,j,1,k)
                                    - T*_S(t,j,0,k)
                                    - (T-1.0)*fac1*_S(t-1,j,0,k)
                                    - (T-2.0)*fac2*_S(t-2,j,0,k)
                                    - (T-3.0)*fac3*_S(t-3,j,0,k);
							}
                        }

						std::for_each(std::execution::par_unseq, grid_idxs.begin(), grid_idxs.end(), [&](uint k){
							_S(t,1,1,k) += 
								recrelS_4(_S(t,1,0), _S(t-1,1,0), _S(t-2,1,0), _S(t-3,1,0), k, p0qq, p1qq, p2qq, p3qq) +
								recrelS_4(_S(t,0,0), _S(t-1,0,0), _S(t-2,0,0), _S(t-3,0,0), k, p0qg, p1qg, p2qg, p3qg);
							_S(t,0,1,k) += 
								recrelS_4(_S(t,1,0), _S(t-1,1,0), _S(t-2,1,0), _S(t-3,1,0), k, p0gq, p1gq, p2gq, p3gq) +
								recrelS_4(_S(t,0,0), _S(t-1,0,0), _S(t-2,0,0), _S(t-3,0,0), k, p0gg, p1gg, p2gg, p3gg);
						});
                    }

                    for (uint t=0; t<=_trunc_idx; ++t) {
                        for (uint j=0; j<=1; ++j) {
                            for (uint k=0; k<_grid.size()-1; k++)
                                arr.get()[j*31][k] += _S(t,j,1,k) * std::pow(_alpha1, t) * std::pow(L1, n)/factorial(n);
							std::ranges::copy(_S(t,j,1), _S(t,j,0).begin());
                        }
                    }
                }
            } break;
        }
		endLogIterations();
    }
	
	void DGLAPSolver::evolveNonSinglet(std::reference_wrapper<std::vector<ArrayGrid>> arr, double L1, double L2, double L3, double L4)
	{
		switch (_order) {
			case 0: { // LO
				log(LOG_DEBUG, "NonSingletLO", "Performing LO non-singlet evolution threaded.");

				std::vector<uint> idxs{};
				for (uint j=13; j<=12+_nf; ++j)
					idxs.emplace_back(j);
				for (uint j=32; j<=30+_nf; ++j)
					idxs.emplace_back(j);
				registerThreadLogs(idxs);

                for (uint j : idxs)
					std::ranges::copy(_A(j,0), arr.get()[j].begin());

				std::vector<std::thread> threads{};
				for (uint j=13; j<=12+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_LO, this, arr, j, L1);	
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_LO, this, arr, j, L1);

			    for (std::thread & t : threads)
					t.join();

				log(LOG_DEBUG, "NonSingletLO", "Finished performing threaded LO non-singlet evolution.");
			} break;
			case 1: { // NLO
				log(LOG_DEBUG, "NonSingletNLO", "Performing NLO non-singlet evolution threaded.");

				std::vector<uint> idxs{};
				for (uint j=13; j<=12+_nf; ++j)
					idxs.emplace_back(j);
				for (uint j=32; j<=30+_nf; ++j)
					idxs.emplace_back(j);
				registerThreadLogs(idxs);

                for (uint j : idxs)
					std::ranges::copy(_B(j,0,0), arr.get()[j].begin());

				std::vector<std::thread> threads{};
                std::array<double, 2> L{L1, L2};
				for (uint j=13; j<=12+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NLO, this, arr, j, ExprName::P1nsm, L);
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NLO, this, arr, j, ExprName::P1nsp, L);

				for (std::thread& t : threads)
					t.join();

				log(LOG_DEBUG, "NonSingletNLO", "Finished performing threaded NLO non-singlet evolution.");
			} break;
			case 2: { // NNLO
				log(LOG_DEBUG, "NonSingletNNLO", "Performing NNLO non-singlet evolution threaded.");

				std::vector<uint> idxs{};
				for (uint j=26; j<=24+_nf; ++j)
					idxs.emplace_back(j);
				for (uint j=32; j<=30+_nf; ++j)
					idxs.emplace_back(j);
				idxs.emplace_back(25);
				registerThreadLogs(idxs);
				
                for (uint j : idxs)
					std::ranges::copy(_C(j,0,0,0), arr.get()[j].begin());

				std::vector<std::thread> threads{};

                std::array<double, 3> L{L1, L2, L3};
				std::array<ExprName, 2> nsm{ExprName::P1nsm, ExprName::P2nsm};
				std::array<ExprName, 2> nsp{ExprName::P1nsp, ExprName::P2nsp};
				std::array<ExprName, 2> nsv{ExprName::P1nsm, ExprName::P2nsv};

				for (uint j=26; j<=24+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NNLO, this, arr, j, nsm, L);
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NNLO, this, arr, j, nsp, L);
				threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NNLO, this, arr, 25, nsv, L);
				
				for (std::thread & t : threads)
					t.join();

				log(LOG_DEBUG, "NonSingletNNLO", "Finished performing threaded NNLO non-singlet evolution.");
			} break;
			case 3: { // N3LO
				log(LOG_DEBUG, "NonSingletN3LO", "Performing N3LO non-singlet evolution threaded.");

				std::vector<uint> idxs{};
				for (uint j=26; j<=24+_nf; ++j)
					idxs.emplace_back(j);
				for (uint j=32; j<=30+_nf; ++j)
					idxs.emplace_back(j);
				idxs.emplace_back(25);
				registerThreadLogs(idxs);

				for (uint j : idxs)
					std::ranges::copy(_D(j,0,0,0,0), arr.get()[j].begin());
				
				std::vector<std::thread> threads{};

                std::array<double, 4> L{L1, L2, L3, L4};
				std::array<ExprName, 3> nsm{ExprName::P1nsm, ExprName::P2nsm, ExprName::P3nsm};
				std::array<ExprName, 3> nsp{ExprName::P1nsp, ExprName::P2nsp, ExprName::P3nsp};
				std::array<ExprName, 3> nsv{ExprName::P1nsm, ExprName::P2nsv, ExprName::P3nsv};

				for (uint j=26; j<=24+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_N3LO, this, arr, j, nsm, L);
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_N3LO, this, arr, j, nsp, L);
				threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_N3LO, this, arr, 25, nsv, L);
				
				for (std::thread& t : threads)
					t.join();

				log(LOG_DEBUG, "NonSingletN3LO", "Finished performing threaded N3LO non-singlet evolution.");
			} break;
		}
	}


    void DGLAPSolver::_mt_EvolveDistribution_NS_LO (
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
        uint j, double L1)
    {
		auto& p0ns = getExpression(ExprName::P0ns);
		
        for (uint n=1; n<_iterations; n++) {
			logThreadIterations(j, n, _iterations-1, "NonSingletLO");
			for (uint k=0; k<_grid.size()-1; k++) {
				_A(j,1,k) = recrelLO(_A(j,0), k, p0ns);
                arr.get()[j][k] += _A(j,1,k)*std::pow(L1, n)/factorial(n);
            }
			std::ranges::copy(_A(j,1), _A(j,0).begin());
		}
    }
    void DGLAPSolver::_mt_EvolveDistribution_NS_NLO(
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
        uint j, ExprName P1, std::array<double, 2> const&  L)
    {
		auto& p0ns = getExpression(ExprName::P0ns);
		auto& p1   = getExpression(P1);
		
        double L1 = L[0];
        double L2 = L[1];
		
        for (uint s=1; s<_iterations; s++) {
            logThreadIterations(j, s, _iterations-1, "NonSingletNLO");
			for (uint n=1; n<=s; n++) {
				for (uint k=0; k<_grid.size()-1;k++) {
					_B(j,1,n,k) = recrelNLO_1(_B(j,0,n-1), k, p0ns);
					arr.get()[j][k] += _B(j,1,n,k)*std::pow(L1,n)*std::pow(L2,s-n)/factorial(n)/factorial(s-n);
				}
			}
			
            for (uint k=0; k<_grid.size()-1;k++) {
                uint n = 0;
                double res = recrelNLO_2(_B(j,0,0), k, p0ns, p1);
                _B(j,1,0,k) = -_B(j,1,1,k) + res;
                arr.get()[j][k] += _B(j,1,0,k)
                    *std::pow(L1,n)*std::pow(L2,s-n)
                    /factorial(n)/factorial(s-n);
            }
			
            for (uint n=0; n<=s; ++n)
				std::ranges::copy(_B(j,1,n), _B(j,0,n).begin());
        }
    }
    void DGLAPSolver::_mt_EvolveDistribution_NS_NNLO(
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
        uint j, std::array<ExprName, 2> const& P, std::array<double, 3> const& L)
    {
		auto& p0ns = getExpression(ExprName::P0ns);
		auto& p1   = getExpression(P[0]);
		auto& p2   = getExpression(P[1]);
		
        double L1 = L[0];
        double L2 = L[1];
        double L3 = L[2];

        for (uint s=1; s<_iterations; s++) {
			logThreadIterations(j, s, _iterations-1, "NonSingletNNLO");
			// recrel #1:
			for (uint t=1; t<=s; t++) {
				for (uint n=1; n<=t; n++) {
					for (uint k=0; k<_grid.size()-1; k++) {
						_C(j,1,t,n,k) = recrelNNLO_1(_C(j,0,t-1,n-1), k, p0ns);

						double orig = _C(j,1,t,n,k);
						double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
						double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
						double res = orig*powers/factorials;

						arr.get()[j][k] += res;
					}
				}
			}

			// recrel #2:
            for (uint k=0; k<_grid.size()-1; k++) {
				double fac1 = -0.5*_C(j,1,s,1,k);
				double fac2 = recrelNNLO_2(_C(j,0,s-1,0), k, p0ns, p1, p2);
				_C(j,1,s,0,k) = fac1 + fac2;

				uint n = 0;
				uint t = s;
				double orig = _C(j,1,s,0,k);
				double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
				double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
				double res = orig*powers/factorials;

				arr.get()[j][k] += res;
			}

			// recrel #3:
			for (int t=s-1; t>=0; t--) {
				for (uint k=0; k<_grid.size()-1; k++) {
					double fac1 = -2.0*_alpha_s.beta1()*(_C(j,1,t+1,0,k) + _C(j,1,t+1,1,k));
					double fac2 = recrelNNLO_3(_C(j,0,t,0), k, p0ns, p1);
					_C(j,1,t,0,k) = fac1 + fac2;

					uint n = 0;
					double orig = _C(j,1,t,0,k);
					double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
					double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
					double res = orig*powers/factorials;

					arr.get()[j][k] += res;
				}
            }

			
            for (uint t=0; t<=s; ++t) {
                for (uint n=0; n<=t; ++n)
					std::ranges::copy(_C(j,1,t,n), _C(j,0,t,n).begin());
            }
        }
    }
    void DGLAPSolver::_mt_EvolveDistribution_NS_N3LO(
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::array<ExprName, 3> const& P, std::array<double, 4> const& L)
    {
		auto& p0ns = getExpression(ExprName::P0ns);
		auto& p1   = getExpression(P[0]);
		auto& p2   = getExpression(P[1]);
		auto& p3   = getExpression(P[2]);
		
		double L1 = L[0];
        double L2 = L[1];
        double L3 = L[2];
        double L4 = L[3];

        // some shorthand
        double r1 = _r1[_nf];
        double b = _b[_nf];
        double c = _c[_nf];
        double gamma = (r1*r1 + r1*b + c)*_alpha_s.beta3();

        for (uint s=1; s<_iterations; s++) {
			logThreadIterations(j, s, _iterations-1, "NonSingletN3LO");

			// recrel #1:
			for (uint t=1; t<=s; t++) {
				for (uint m=1; m<=t; m++) {
					for (uint n=1; n<=m; n++) {
						for (uint k=0; k<_grid.size()-1; k++) {
							_D(j,1,t,m,n,k) = recrelN3LO_1(_D(j,0,t-1,m-1,n-1), k, p0ns);

							double orig = _D(j,1,t,m,n,k);
							double powers =
								std::pow(L1,n)
								*std::pow(L2,(m-n))
								*std::pow(L3,(t-m))
								*std::pow(L4,(s-t));
							double factorials =
								factorial(n)
								*factorial(m-n)
								*factorial(t-m)
								*factorial(s-t);
							double res = orig*powers/factorials;
                            
							arr.get()[j][k] += res;
						}
					}
				}
			}

			// recrel #2:
            for (uint k=0; k<_grid.size()-1; k++) {
				double fac1 = (
					0.5*(16.0*PI_2*_alpha_s.beta1() + 4*PI*r1*_alpha_s.beta2() - (c + b*r1)*_alpha_s.beta3())
				) * _D(j,1,s,s,1,k);
				double fac2 = recrelN3LO_2(_D(j,0,s-1,s-1,0), k, p0ns, p1, p2, p3);
				_D(j,1,s,s,0,k) = (fac1 + fac2)/gamma;

				uint t = s;
				uint m = s;
				uint n = 0;
				double orig = _D(j,1,s,s,0,k);
				double powers =
					std::pow(L1,n)
					*std::pow(L2,(m-n))
					*std::pow(L3,(t-m))
					*std::pow(L4,(s-t));
				double factorials =
					factorial(n)
					*factorial(m-n)
					*factorial(t-m)
					*factorial(s-t);
				double res = orig*powers/factorials;
                    
				arr.get()[j][k] += res;
			}

			// recrel #3:
			for (int m=s-1; m>=0; m--) {
				for (uint k=0; k<_grid.size()-1; k++) {
                    double fac1 = -(
                        16.0*PI_2*_alpha_s.beta1() + 4.0*PI*r1*_alpha_s.beta2() + r1*r1*_alpha_s.beta3()
                    ) * _D(j,1,s,m+1,1,k);
                    double fac2 = recrelN3LO_3(_D(j,0,s-1,m,0), k, p0ns, p1, p2, p3);
                    _D(j,1,s,m,0,k) = (fac1 + fac2)/gamma;

                    uint t = s;
                    uint n = 0;
                    double orig = _D(j,1,s,m,0,k);
                    double powers =
                        std::pow(L1,n)
                        *std::pow(L2,(m-n))
                        *std::pow(L3,(t-m))
                        *std::pow(L4,(s-t));
                    double factorials =
                        factorial(n)
                        *factorial(m-n)
                        *factorial(t-m)
                        *factorial(s-t);
                    double res = orig*powers/factorials;
                    
                    arr.get()[j][k] += res;
                }
			}

			// recrel #4:
			for (int t=s-1; t>=0; t--) {
				for (int m=t; m>=0; m--) {
					for (uint k=0; k<_grid.size()-1; k++) {
						double fac1a = -2.0*b*gamma;
						double fac1b = 32*PI_2*(b+r1)*_alpha_s.beta1() - 8*PI*c*_alpha_s.beta2() - 2*c*r1*_alpha_s.beta3();
						double fac1 = fac1a*_D(j,1,t+1,m+1,0,k) + fac1b*_D(j,1,t+1,m+1,1,k);
						double fac2 = recrelN3LO_4(_D(j,0,t,m,0), k, p0ns, p1, p2, p3);
						_D(j,1,t,m,0,k) = (fac1 + fac2)/gamma;

						uint n = 0;
						double orig = _D(j,1,t,m,0,k);
						double powers =
							std::pow(L1,n)
							*std::pow(L2,(m-n))
							*std::pow(L3,(t-m))
							*std::pow(L4,(s-t));
						double factorials =
							factorial(n)
							*factorial(m-n)
							*factorial(t-m)
							*factorial(s-t);
						double res = orig*powers/factorials;
                        
						arr.get()[j][k] += res;
					}
				}
			}

            for (uint t=0; t<=s; ++t) {
                for (uint m=0; m<=s; ++m) {
                    for (uint n=0; n<=s; ++n)
						std::ranges::copy(_D(j,1,t,m,n), _D(j,0,t,m,n).begin());
                }
            }
        }
    }

	void DGLAPSolver::evolveNonSingletTrunc(std::reference_wrapper<std::vector<ArrayGrid>> arr, double L1)
	{
		log(LOG_DEBUG, "DGLAP", "Using the truncated ansatz for the non-singlet sector");
        switch (_order) {
            case 0: {
				log(LOG_DEBUG, "NonSingletLO (trunc)", "Performing LO non-singlet evolution threaded.");
				
				std::vector<uint> idxs{};
				for (uint j=13; j<=12+_nf; ++j)
					idxs.emplace_back(j);
				for (uint j=32; j<=30+_nf; ++j)
					idxs.emplace_back(j);
				registerThreadLogs(idxs);

				for (uint j : idxs)
					std::ranges::copy(_S_NS(0,j,0), arr.get()[j].begin());

				std::vector<std::thread> threads{};
				for (uint j=13; j<=12+_nf; ++j)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NST_LO, this, arr, j, L1);
				for (uint j=32; j<=30+_nf; ++j)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NST_LO, this, arr, j, L1);

				for (std::thread& t : threads)
					t.join();
				
				log(LOG_DEBUG, "NonSingletLO (trunc)", "Finished performing LO non-singlet evolution threaded.");
            }; break;
			case 1: {
				log(LOG_DEBUG, "NonSingletNLO (trunc)", "Performing NLO non-singlet evolution threaded.");

				std::vector<uint> idxs{};
				for (uint j=13; j<=12+_nf; ++j)
					idxs.emplace_back(j);
				for (uint j=32; j<=30+_nf; ++j)
					idxs.emplace_back(j);
				registerThreadLogs(idxs);

				for (uint j : idxs)
					std::ranges::copy(_S_NS(0,j,0), arr.get()[j].begin());

				std::vector<std::thread> threads{};
				for (uint j=13; j<=12+_nf; ++j)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NST_NLO, this, arr, j, ExprName::P1nsm, L1);
				for (uint j=32; j<=30+_nf; ++j)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NST_NLO, this, arr, j, ExprName::P1nsp, L1);

				for (std::thread & t : threads)
					t.join();
				
				log(LOG_DEBUG, "NonSingletNLO (trunc)", "Finished performing NLO non-singlet evolution threaded.");
			}; break;
			case 2: {
				log(LOG_DEBUG, "NonSingletNNLO (trunc)", "Performing NNLO non-singlet evolution threaded.");

				std::vector<uint> idxs{};
				for (uint j=26; j<=24+_nf; j++)
					idxs.emplace_back(j);
				for (uint j=32; j<=30+_nf; ++j)
					idxs.emplace_back(j);
				idxs.emplace_back(25);
				registerThreadLogs(idxs);

				for (uint j : idxs)
					std::ranges::copy(_S_NS(0,j,0), arr.get()[j].begin());
				
				std::vector<std::thread> threads{};
				for (uint j=26; j<=24+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NST_NNLO, this, arr, j, ExprName::P1nsm, ExprName::P2nsm, L1);
				for (uint j=32; j<=30+_nf; ++j)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NST_NNLO, this, arr, j, ExprName::P1nsp, ExprName::P2nsp, L1);
				threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NST_NNLO, this, arr, 25, ExprName::P1nsm, ExprName::P2nsv, L1);
				
				for (std::thread & t : threads)
					t.join();
				
				log(LOG_DEBUG, "NonSingletNNLO (trunc)", "Finished performing NNLO non-singlet evolution threaded.");
			}; break;
			case 3: {
				log(LOG_DEBUG, "NonSingletN3LO (trunc)", "Performing N3LO non-singlet evolution threaded.");

				std::vector<uint> idxs{};
				for (uint j=26; j<=24+_nf; j++)
					idxs.emplace_back(j);
				for (uint j=32; j<=30+_nf; ++j)
					idxs.emplace_back(j);
				idxs.emplace_back(25);
				registerThreadLogs(idxs);

				for (uint j : idxs)
					std::ranges::copy(_S_NS(0,j,0), arr.get()[j].begin());
				
				std::vector<std::thread> threads{};
				for (uint j=26; j<=24+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NST_N3LO, this, arr, j, ExprName::P1nsm, ExprName::P2nsm, ExprName::P3nsm, L1);
				for (uint j=32; j<=30+_nf; ++j)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NST_N3LO, this, arr, j, ExprName::P1nsp, ExprName::P2nsp, ExprName::P3nsp, L1);
						threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NST_N3LO, this, arr, 25, ExprName::P1nsm, ExprName::P2nsv, ExprName::P3nsv, L1);

				for (std::thread & t : threads)
					t.join();
				
				log(LOG_DEBUG, "NonSingletN3LO (trunc)", "Finished performing N3LO non-singlet evolution threaded.");
			}; break;
        }
	}

	void DGLAPSolver::_mt_EvolveDistribution_NST_LO(std::reference_wrapper<std::vector<ArrayGrid>> arr, uint j, double L1)
	{
		auto& p0ns = getExpression(ExprName::P0ns);

		for (uint n=1; n<_iterations; n++) {
			logThreadIterations(j, n, _iterations-1, "NonSingletLO (trunc)");
			double fac = std::pow(L1, n)/factorial(n);

			for (uint k=0; k<_grid.size()-1; k++) {
				_S_NS(0,j,1,k) = recrelS_1(_S_NS(0,j,0), k, p0ns);

				arr.get()[j][k] += _S_NS(0,j,1,k)*fac;
			}

			std::ranges::copy(_S_NS(0,j,1), _S_NS(0,j,0).begin());
		}
	}
	void DGLAPSolver::_mt_EvolveDistribution_NST_NLO(
		std::reference_wrapper<std::vector<ArrayGrid>> arr, uint j,
		ExprName _p1, double L1)
	{

		auto& p0ns = getExpression(ExprName::P0ns);
		auto& p1 =   getExpression(_p1);

		for (uint n=1; n<_iterations; n++) {
			logThreadIterations(j, n, _iterations-1, "NonSingletNLO (trunc)");
			double fac = std::pow(L1, n)/factorial(n);

			// LO piece
			for (uint k=0; k<_grid.size()-1; k++)
				_S_NS(0,j,1,k) = recrelS_1(_S_NS(0,j,0), k, p0ns);

			// NLO
			for (uint k=0; k<_grid.size()-1; k++) {
				_S_NS(1,j,1,k) = -_S_NS(0,j,1,k) * _alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()) - _S_NS(1,j,0,k);
				_S_NS(1,j,1,k) += recrelS_2(_S_NS(1,j,0), _S_NS(0,j,0), k, p0ns, p1);
			}

			// truncation terms
			for (uint t=2; t<=_trunc_idx; ++t) {
				double T = static_cast<double>(t);
				for (uint k=0; k<_grid.size()-1; k++) {
					_S_NS(t,j,1,k) =
						- (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S_NS(t-1,j,1,k)
						- T*_S_NS(t,j,0,k)
						- (T-1.0)*(_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S_NS(t-1,j,0,k);
					_S_NS(t,j,1,k) += recrelS_2(_S_NS(t,j,0), _S_NS(t-1,j,0), k, p0ns, p1);
				}
			}

			// resum
			for (uint t=0; t<=_trunc_idx; ++t) {
				double a = std::pow(_alpha1, t);
				for (uint k=0; k<_grid.size()-1; k++)
					arr.get()[j][k] += _S_NS(t,j,1,k)*a*fac;
			}

			// setup for next iteration
			for (uint t=0; t<=_trunc_idx; ++t)
				std::ranges::copy(_S_NS(t,j,1), _S_NS(t,j,0).begin());
		}
	}
	void DGLAPSolver::_mt_EvolveDistribution_NST_NNLO(std::reference_wrapper<std::vector<ArrayGrid>> arr, uint j,
		ExprName _p1, ExprName _p2, double L1)
	{
		auto& p0ns = getExpression(ExprName::P0ns);
		auto& p1 = getExpression(_p1);
		auto& p2 = getExpression(_p2);

		for (uint n=1; n<_iterations; n++) {
			logThreadIterations(j, n, _iterations-1, "NonSingletNNLO (trunc)");
			double fac = std::pow(L1, n)/factorial(n);

			// LO piece
			for (uint k=0; k<_grid.size()-1; k++)
				_S_NS(0,j,1,k) = recrelS_1(_S_NS(0,j,0), k, p0ns);

			// NLO
			for (uint k=0; k<_grid.size()-1; k++) {
				_S_NS(1,j,1,k) = -_S_NS(0,j,1,k) * _alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()) - _S_NS(1,j,0,k);
				_S_NS(1,j,1,k) += recrelS_2(_S_NS(1,j,0), _S_NS(0,j,0), k, p0ns, p1);
			}

			// NNLO
			for (uint k=0; k<_grid.size()-1; k++) {
				_S_NS(2,j,1,k) =
					- (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S_NS(1,j,1,k)
					- (_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S_NS(0,j,1,k)
					- 2.0*_S_NS(2,j,0,k)
					- (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S_NS(1,j,0,k);
				_S_NS(2,j,1,k) += recrelS_3(_S_NS(2,j,0), _S_NS(1,j,0), _S_NS(0,j,0), k, p0ns, p1, p2);
			}

			for (uint t=3; t<=_trunc_idx; ++t) {
				double T = static_cast<double>(t);
				for (uint k=0; k<_grid.size()-1; k++) {
					_S_NS(t,j,1,k) =
						- (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S_NS(t-1,j,1,k)
						- (_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S_NS(t-2,j,1,k)
						- T*_S_NS(t,j,0,k)
						- (T-1.0)*(_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S_NS(t-1,j,0,k)
						- (T-2.0)*(_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S_NS(t-2,j,0,k);
					_S_NS(t,j,1,k) += recrelS_3(_S_NS(t,j,0), _S_NS(t-1,j,0), _S_NS(t-2,j,0), k, p0ns, p1, p2);
				}
			}

			// resum
			for (uint t=0; t<=_trunc_idx; ++t) {
				double a = std::pow(_alpha1, t);
				for (uint k=0; k<_grid.size()-1; k++)
					arr.get()[j][k] += _S_NS(t,j,1,k)*a*fac;
			}

			// setup for next iteration
			for (uint t=0; t<=_trunc_idx; ++t)
				std::ranges::copy(_S_NS(t,j,1), _S_NS(t,j,0).begin());
		}
	}
	void DGLAPSolver::_mt_EvolveDistribution_NST_N3LO(
		std::reference_wrapper<std::vector<ArrayGrid>> arr, uint j,
		ExprName _p1, ExprName _p2, ExprName _p3, double L1)
	{
		auto& p0ns = getExpression(ExprName::P0ns);
		auto& p1 =   getExpression(_p1);
		auto& p2 =   getExpression(_p2);
		auto& p3 =   getExpression(_p3);

		double beta0 = _alpha_s.beta0();
		double beta1 = _alpha_s.beta1();
		double beta2 = _alpha_s.beta2();
		double beta3 = _alpha_s.beta3();

		for (uint n=1; n<_iterations; n++) {
			logThreadIterations(j, n, _iterations-1, "NonSingletN3LO (trunc)");
			double fac = std::pow(L1, n)/factorial(n);

			// LO piece
			for (uint k=0; k<_grid.size()-1; k++)
				_S_NS(0,j,1,k) = recrelS_1(_S_NS(0,j,0), k, p0ns);

			// NLO
			for (uint k=0; k<_grid.size()-1; k++) {
				_S_NS(1,j,1,k) = -_S_NS(0,j,1,k) * beta1/(4.0*PI*beta0) - _S_NS(1,j,0,k);
				_S_NS(1,j,1,k) += recrelS_2(_S_NS(1,j,0), _S_NS(0,j,0), k, p0ns, p1);
			}

			// NNLO
			for (uint k=0; k<_grid.size()-1; k++) {
				_S_NS(2,j,1,k) =
					- (beta1/(4.0*PI*beta0))*_S_NS(1,j,1,k)
					- (beta2/(16.0*PI_2*beta0))*_S_NS(0,j,1,k)
					- 2.0*_S_NS(2,j,0,k)
					- (beta1/(4.0*PI*beta0))*_S_NS(1,j,0,k);
				_S_NS(2,j,1,k) += recrelS_3(_S_NS(2,j,0), _S_NS(1,j,0), _S_NS(0,j,0), k, p0ns, p1, p2);
			}

			// N3LO
			for (uint k=0; k<_grid.size()-1; k++) {
				_S_NS(3,j,1,k) =
					- (beta1/(4.0*PI*beta0))*_S_NS(2,j,1,k)
					- (beta2/(16.0*PI_2*beta0))*_S_NS(1,j,1,k)
					- (beta3/(64.0*PI_3*beta0))*_S_NS(0,j,1,k)
					- 3.0*_S_NS(3,j,0,k)
					- 2.0*(beta1/(4.0*PI*beta0))*_S_NS(2,j,0,k)
					- (beta2/(16.0*PI_2*beta0))*_S_NS(1,j,0,k);
				_S_NS(3,j,1,k) += recrelS_4(
					_S_NS(3,j,0), _S_NS(2,j,0), _S_NS(1,j,0), _S_NS(0,j,0),
					k,
					p0ns, p1, p2, p3);
			}

			for (uint t=4; t<=_trunc_idx; ++t) {
				double T = static_cast<double>(t);
				for (uint k=0; k<_grid.size()-1; k++) {
					_S_NS(t,j,1,k) =
						- (beta1/(4.0*PI*beta0))*_S_NS(t-1,j,1,k)
						- (beta2/(16.0*PI_2*beta0))*_S_NS(t-2,j,1,k)
						- (beta3/(64.0*PI_3*beta0))*_S_NS(t-3,j,1,k)
						- T*_S_NS(t,j,0,k)
						- (T-1.0)*(beta1/(4.0*PI*beta0))*_S_NS(t-1,j,0,k)
						- (T-2.0)*(beta2/(16.0*PI_2*beta0))*_S_NS(t-2,j,0,k)
						- (T-3.0)*(beta3/(64.0*PI_3*beta0))*_S_NS(t-3,j,0,k);
					_S_NS(t,j,1,k) += recrelS_4(
						_S_NS(t,j,0), _S_NS(t-1,j,0), _S_NS(t-2,j,0), _S_NS(t-3,j,0),
						k,
						p0ns, p1, p2, p3);
				}
			}
					
			// resum
			for (uint t=0; t<=_trunc_idx; ++t) {
				double a = std::pow(_alpha1, t);
				for (uint k=0; k<_grid.size()-1; k++)
					arr.get()[j][k] += _S_NS(t,j,1,k)*a*fac;
			}

			// setup for next iteration
			for (uint t=0; t<=_trunc_idx; ++t)
				std::ranges::copy(_S_NS(t,j,1), _S_NS(t,j,0).begin());
		}
	}	
} // namespace Candia2
