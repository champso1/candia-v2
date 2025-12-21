#include "Candia-v2/Candia.hpp"
#include "Candia-v2/FuncArrayGrid.hpp"
#include "Candia-v2/Math.hpp"

#include <cmath>
#include <functional>
#include <print>

namespace Candia2
{
    void DGLAPSolver::evolveSinglet(
        std::reference_wrapper<std::vector<ArrayGrid>> arr, double L1
    ) {
        for (uint j=0; j<=1; ++j)
            arr.get()[j*31] = _S2[0][j][0];

        switch (_order)
        {
            case 0:
            {
                for (uint n=1; n<_iterations; n++)
                {
                    std::println("[DGLAP] LO Singlet Iteration {}", n);
                    
                    for (uint k=0; k<_grid.size()-1; k++)
                    {
                        _S2[0][1][1][k] =
							recrelS_1(_S2[0][1][0], k, getExpression("P0qq")) +
							recrelS_1(_S2[0][0][0], k, getExpression("P0qg"));
                        _S2[0][0][1][k] =
							recrelS_1(_S2[0][1][0], k, getExpression("P0gq")) +
							recrelS_1(_S2[0][0][0], k, getExpression("P0gg"));

						for (uint j=0; j<=1; j++)
							arr.get()[j*31][k] += _S2[0][j][1][k] * std::pow(L1, n)/factorial(n);
                    }

                    for (uint j=0; j<=1; ++j)
                        _S2[0][j][0] = _S2[0][j][1];
                }
            } break;
            case 1:
            {
                for (uint n=1; n<_iterations; n++)
                {
                    std::println("[DGLAP] NLO Singlet Iteration {}", n);
                    
                    for (uint k=0; k<_grid.size()-1; k++)
                    {
						// LO piece
                        _S2[0][1][1][k] =
							recrelS_1(_S2[0][1][0], k, getExpression("P0qq")) +
							recrelS_1(_S2[0][0][0], k, getExpression("P0qg"));
                        _S2[0][0][1][k] =
							recrelS_1(_S2[0][1][0], k, getExpression("P0gq")) +
							recrelS_1(_S2[0][0][0], k, getExpression("P0gg"));

						// NLO convolution piece
						_S2[1][1][1][k] +=
							recrelS_2(_S2[0][1][0],_S2[1][1][0],k,getExpression("P0qq"),getExpression("P1qq")) +
							recrelS_2(_S2[0][0][0],_S2[1][0][0],k,getExpression("P0qg"),getExpression("P1qg"));
					    _S2[1][0][1][k] +=
							recrelS_2(_S2[0][1][0],_S2[1][1][0],k,getExpression("P0gq"),getExpression("P1gq"))+
							recrelS_2(_S2[0][0][0],_S2[1][0][0],k,getExpression("P0gg"),getExpression("P1gg"));

						// other NLO piece
						for (uint j=0; j<=1; j++)
                            _S2[1][j][1][k] = -_S2[0][j][1][k] * _alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()) - _S2[1][j][0][k];
                    }

                    // NLO truncation terms
                    for (uint t=2; t<=_trunc_idx; ++t)
                    {
                        double T = static_cast<double>(t);

                        // non-convolution piece:
                        for (uint k=0; k<_grid.size()-1; k++)
                        {
                            for (uint j=0; j<=1; j++)
                                _S2[t][j][1][k] =
                                    - (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S2[t-1][j][1][k]
                                    - T*_S2[t][j][0][k]
                                    - (T-1.0)*(_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S2[t-1][j][0][k];
                        }

                        // convolution piece
                        for (uint k=0; k<_grid.size()-1; k++)
                        {
                            _S2[t][1][1][k] +=
								recrelS_2(_S2[t-1][1][0],_S2[t][1][0],k,getExpression("P0qq"),getExpression("P1qq"))+
								recrelS_2(_S2[t-1][0][0],_S2[t][0][0],k,getExpression("P0qg"),getExpression("P1qg"));
                            _S2[t][0][1][k] +=
								recrelS_2(_S2[t-1][1][0],_S2[t][1][0],k,getExpression("P0gq"),getExpression("P1gq"))+
								recrelS_2(_S2[t-1][0][0],_S2[t][0][0],k,getExpression("P0gg"),getExpression("P1gg"));
                        }
                    }

                    for (uint t=0; t<=_trunc_idx; ++t)
                    {
                        for (uint j=0; j<=1; ++j)
                        {
                            for (uint k=0; k<_grid.size()-1; k++)
                                arr.get()[j*31][k] += _S2[t][j][1][k] * std::pow(_alpha1, t) * std::pow(L1, n)/factorial(n);        
                        }
                    }
                
                    for (uint t=0; t<=_trunc_idx; ++t)
                    {
                        for (uint j=0; j<=1; ++j)
                            _S2[t][j][0] = _S2[t][j][1];
                    }
                }
            } break;
            case 2:
            {
                for (uint n=1; n<_iterations; n++)
                {
                    std::println("[DGLAP] NNLO Singlet Iteration {}", n);

                    for (uint k=0; k<_grid.size()-1; k++)
                    {
						// LO piece
                        _S2[0][1][1][k] =
							recrelS_1(_S2[0][1][0], k, getExpression("P0qq")) +
							recrelS_1(_S2[0][0][0], k, getExpression("P0qg"));
                        _S2[0][0][1][k] =
							recrelS_1(_S2[0][1][0], k, getExpression("P0gq")) +
							recrelS_1(_S2[0][0][0], k, getExpression("P0gg"));

						// NLO convolution piece
						_S2[1][1][1][k] +=
							recrelS_2(_S2[0][1][0],_S2[1][1][0],k,getExpression("P0qq"),getExpression("P1qq"))+
							recrelS_2(_S2[0][0][0],_S2[1][0][0],k,getExpression("P0qg"),getExpression("P1qg"));
					    _S2[1][0][1][k] +=
							recrelS_2(_S2[0][1][0],_S2[1][1][0],k,getExpression("P0gq"),getExpression("P1gq"))+
							recrelS_2(_S2[0][0][0],_S2[1][0][0],k,getExpression("P0gg"),getExpression("P1gg"));

						// NNLO convolution piece
						_S2[2][1][1][k] += 
                            recrelS_3(
                                _S2[0][1][0], _S2[1][1][0], _S2[2][1][0], k, 
                                getExpression("P0qq"), getExpression("P1qq"), getExpression("P2qq")) +
                            recrelS_3(
                                _S2[0][0][0], _S2[1][0][0], _S2[2][0][0], k,
                                getExpression("P0qg"), getExpression("P1qg"), getExpression("P2qg"));

                        _S2[2][0][1][k] += 
                            recrelS_3(_S2[0][1][0], _S2[1][1][0], _S2[2][1][0], k,
                                getExpression("P0gq"), getExpression("P1gq"), getExpression("P2gq")) +
                            recrelS_3(_S2[0][0][0], _S2[1][0][0], _S2[2][0][0], k,
                                getExpression("P0gg"), getExpression("P1gg"), getExpression("P2gg"));

						// other NLO piece
						for (uint j=0; j<=1; j++)
                            _S2[1][j][1][k] = -_S2[0][j][1][k] * _alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()) - _S2[1][j][0][k];

						// other NNLO piece
						for (uint j=0; j<=1; j++)
						{
                            _S2[2][j][1][k] =
                                - (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S2[1][j][1][k]
                                - (_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S2[0][j][1][k]
                                - 2.0*_S2[2][j][0][k]
                                - (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S2[1][j][0][k];
						}

						
                    }
					
                    for (uint t=3; t<=_trunc_idx; ++t)
                    {
                        double T = static_cast<double>(t);

                        // non-convolution piece:
                        for (uint k=0; k<_grid.size()-1; k++)
                        {
                            for (uint j=0; j<=1; j++)
                                _S2[t][j][1][k] =
                                    - (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S2[t-1][j][1][k]
                                    - (_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S2[t-2][j][1][k]
                                    - T*_S2[t][j][0][k]
                                    - (T-1.0)*(_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S2[t-1][j][0][k]
                                    - (T-2.0)*(_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S2[t-2][j][0][k];
                        }

                        // convolution piece
                        for (uint k=0; k<_grid.size()-1; k++)
                        {
                            _S2[t][1][1][k] +=
								recrelS_3(
									_S2[t-2][1][0], _S2[t-1][1][0], _S2[t][1][0], k,
									getExpression("P0qq"), getExpression("P1qq"), getExpression("P2qq")) +
                                recrelS_3(
									_S2[t-2][0][0], _S2[t-1][0][0], _S2[t][0][0], k,
                                    getExpression("P0qg"), getExpression("P1qg"), getExpression("P2qg"));

                            _S2[t][0][1][k] += 
                                recrelS_3(
									_S2[t-2][1][0], _S2[t-1][1][0], _S2[t][1][0], k,
                                    getExpression("P0gq"), getExpression("P1gq"), getExpression("P2gq")) +
                                recrelS_3(
									_S2[t-2][0][0], _S2[t-1][0][0], _S2[t][0][0], k,
                                    getExpression("P0gg"), getExpression("P1gg"), getExpression("P2gg"));
                        }
                    }

                    for (uint t=0; t<=_trunc_idx; ++t)
                    {
                        for (uint j=0; j<=1; ++j)
                        {
                            for (uint k=0; k<_grid.size()-1; k++)
                                arr.get()[j*31][k] += _S2[t][j][1][k] * std::pow(_alpha1, t) * std::pow(L1, n)/factorial(n);
                        }
                    }
                
                    for (uint t=0; t<=_trunc_idx; ++t)
                    {
                        for (uint j=0; j<=1; ++j)
                            _S2[t][j][0] = _S2[t][j][1];
                    }
                }
            } break;
            case 3:
            {
                for (uint n=1; n<_iterations; n++)
                {
                    std::println("[DGLAP] N3LO Singlet Iteration {}", n);

                    for (uint k=0; k<_grid.size()-1; k++)
                    {
						// LO piece
                        _S2[0][1][1][k] =
							recrelS_1(_S2[0][1][0], k, getExpression("P0qq")) +
							recrelS_1(_S2[0][0][0], k, getExpression("P0qg"));
                        _S2[0][0][1][k] =
							recrelS_1(_S2[0][1][0], k, getExpression("P0gq")) +
							recrelS_1(_S2[0][0][0], k, getExpression("P0gg"));

						// NLO convolution piece
						_S2[1][1][1][k] +=
							recrelS_2(_S2[0][1][0],_S2[1][1][0],k,getExpression("P0qq"),getExpression("P1qq"))+
							recrelS_2(_S2[0][0][0],_S2[1][0][0],k,getExpression("P0qg"),getExpression("P1qg"));
					    _S2[1][0][1][k] +=
							recrelS_2(_S2[0][1][0],_S2[1][1][0],k,getExpression("P0gq"),getExpression("P1gq"))+
							recrelS_2(_S2[0][0][0],_S2[1][0][0],k,getExpression("P0gg"),getExpression("P1gg"));

						// NNLO convolution piece
						_S2[2][1][1][k] += 
                            recrelS_3(
                                _S2[0][1][0], _S2[1][1][0], _S2[2][1][0], k, 
                                getExpression("P0qq"), getExpression("P1qq"), getExpression("P2qq")) +
                            recrelS_3(
                                _S2[0][0][0], _S2[1][0][0], _S2[2][0][0], k,
                                getExpression("P0qg"), getExpression("P1qg"), getExpression("P2qg"));
                        _S2[2][0][1][k] += 
                            recrelS_3(
								_S2[0][1][0], _S2[1][1][0], _S2[2][1][0], k,
                                getExpression("P0gq"), getExpression("P1gq"), getExpression("P2gq")) +
                            recrelS_3(
								_S2[0][0][0], _S2[1][0][0], _S2[2][0][0], k,
                                getExpression("P0gg"), getExpression("P1gg"), getExpression("P2gg"));

						// N3LO convolution piece
						_S2[3][1][1][k] += 
                            recrelS_4(
								_S2[0][1][0], _S2[1][1][0], _S2[2][1][0], _S2[3][1][0], k,
                                getExpression("P0qq"), getExpression("P1qq"),
								getExpression("P2qq"), getExpression("P3qq")) +
                            recrelS_4(
								_S2[0][0][0], _S2[1][0][0], _S2[2][0][0], _S2[3][0][0], k,
                                getExpression("P0qg"), getExpression("P1qg"),
								getExpression("P2qg"), getExpression("P3qg"));
                        _S2[3][0][1][k] += 
                            recrelS_4(
								_S2[0][1][0], _S2[1][1][0], _S2[2][1][0], _S2[3][1][0], k,
                                getExpression("P0gq"), getExpression("P1gq"),
								getExpression("P2gq"), getExpression("P3gq")) +
                            recrelS_4(
								_S2[0][0][0], _S2[1][0][0], _S2[2][0][0], _S2[3][0][0], k,
                                getExpression("P0gg"), getExpression("P1gg"),
								getExpression("P2gg"), getExpression("P3gg"));

						// other NLO piece
						for (uint j=0; j<=1; j++)
                            _S2[1][j][1][k] = -_S2[0][j][1][k] * _alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()) - _S2[1][j][0][k];

						// other NNLO piece
						for (uint j=0; j<=1; j++)
						{
                            _S2[2][j][1][k] =
                                - (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S2[1][j][1][k]
                                - (_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S2[0][j][1][k]
                                - 2.0*_S2[2][j][0][k]
                                - (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S2[1][j][0][k];
						}

						// other N3LO piece
						for (uint j=0; j<=1; j++)
						{
                            _S2[3][j][1][k] =
                                - (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S2[2][j][1][k]
                                - (_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S2[1][j][1][k]
                                - (_alpha_s.beta3()/(64.0*PI_3*_alpha_s.beta0()))*_S2[0][j][1][k]
                                - 3.0*_S2[3][j][0][k]
                                - 2.0*(_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S2[2][j][0][k]
                                - (_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S2[1][j][0][k];
						}
                    }

                    for (uint t=4; t<=_trunc_idx; ++t)
                    {
                        double T = static_cast<double>(t);

                        // non-convolution piece:
                        for (uint k=0; k<_grid.size()-1; k++)
                        {
                            for (uint j=0; j<=1; j++)
                                _S2[t][j][1][k] =
                                    - (_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S2[t-1][j][1][k]
                                    - (_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S2[t-2][j][1][k]
                                    - (_alpha_s.beta3()/(64.0*PI_3*_alpha_s.beta0()))*_S2[t-3][j][1][k]
                                    - T*_S2[t][j][0][k]
                                    - (T-1.0)*(_alpha_s.beta1()/(4.0*PI*_alpha_s.beta0()))*_S2[t-1][j][0][k]
                                    - (T-2.0)*(_alpha_s.beta2()/(16.0*PI_2*_alpha_s.beta0()))*_S2[t-2][j][0][k]
                                    - (T-3.0)*(_alpha_s.beta3()/(64.0*PI_3*_alpha_s.beta0()))*_S2[t-3][j][0][k];
                        }

                        // convolution piece
                        for (uint k=0; k<_grid.size()-1; k++)
                        {
                            _S2[t][1][1][k] += 
                                recrelS_4(
									_S2[t-3][1][0], _S2[t-2][1][0], _S2[t-1][1][0], _S2[t][1][0],k,
                                    getExpression("P0qq"), getExpression("P1qq"),
									getExpression("P2qq"), getExpression("P3qq")) +
                                recrelS_4(
									_S2[t-3][0][0], _S2[t-2][0][0], _S2[t-1][0][0], _S2[t][0][0],k,
                                    getExpression("P0qg"), getExpression("P1qg"),
									getExpression("P2qg"), getExpression("P3qg"));
                            _S2[t][0][1][k] += 
                                recrelS_4(
									_S2[t-3][1][0], _S2[t-2][1][0], _S2[t-1][1][0], _S2[t][1][0],k,
                                    getExpression("P0gq"), getExpression("P1gq"),
									getExpression("P2gq"), getExpression("P3gq")) +
                                recrelS_4(
									_S2[t-3][0][0], _S2[t-2][0][0], _S2[t-1][0][0], _S2[t][0][0],k,
                                    getExpression("P0gg"), getExpression("P1gg"),
									getExpression("P2gg"), getExpression("P3gg"));
                        }
                    }

                    for (uint t=0; t<=_trunc_idx; ++t)
                    {
                        for (uint j=0; j<=1; ++j)
                        {
                            for (uint k=0; k<_grid.size()-1; k++)
                                arr.get()[j*31][k] += _S2[t][j][1][k] * std::pow(_alpha1, t) * std::pow(L1, n)/factorial(n);        
                        }
                    }
                
                    for (uint t=0; t<=_trunc_idx; ++t)
                    {
                        for (uint j=0; j<=1; ++j)
                            _S2[t][j][0] = _S2[t][j][1];
                    }
                }
            } break;
        }
    }



    void DGLAPSolver::evolveNonSinglet(
        std::reference_wrapper<std::vector<ArrayGrid>> arr,
        double L1, double L2, double L3, double L4) 
    {
        switch (_order)
        {
            case 0:
            {
                for (uint j=13; j<=12+_nf; ++j)
                    arr.get()[j] = _A2[j][0];
                for (uint j=32; j<=30+_nf; ++j)
                    arr.get()[j] = _A2[j][0];

                for (uint n=1; n<_iterations; n++)
                {
                    std::println("[DGLAP] LO Non-Singlet Iteration {}", n);
                    for (uint k=0; k<_grid.size()-1; k++)
                    {
                        for (uint j=13; j<=30+_nf; j++)
                        {
                            _A2[j][1][k] = recrelLO(_A2[j][0], k, getExpression("P0ns"));
                            arr.get()[j][k] += _A2[j][1][k]*std::pow(L1, n)/factorial(n);

                            if (j == (12+_nf))
                                j = 31;
                        }   
                    }

                    for (uint j=13; j<=30+_nf; ++j)
                    {
                        _A2[j][0] = _A2[j][1];
                        if (j == (12+_nf))
                            j = 31;
                    }
                }
            } break;
            case 1:
            {
                for (uint j=13; j<=12+_nf; ++j)
                    arr.get()[j] = _B2[j][0][0];
                for (uint j=32; j<=30+_nf; ++j)
                    arr.get()[j] = _B2[j][0][0];

                for (uint s=1; s<_iterations; s++)
				{
					std::println("[DGLAP] NLO Non-Singlet Iteration {}", s);

					for (uint k=0; k<_grid.size()-1;k++)
					{
						for (uint j=13; j<=12+_nf; j++)
						{
							for (uint n=1; n<=s; n++)
                            {
								_B2[j][1][n][k] = recrelNLO_1(_B2[j][0][n-1], k, getExpression("P0ns"));
                                arr.get()[j][k] += _B2[j][1][n][k]*std::pow(L1,n)*std::pow(L2,s-n)/factorial(n)/factorial(s-n);
                            }

                            uint n = 0;
                            double res = recrelNLO_2(_B2[j][0][0], k, 
                                getExpression("P0ns"), 
                                getExpression("P1nsm"));
							_B2[j][1][0][k] = -_B2[j][1][1][k] + res;
                            arr.get()[j][k] += _B2[j][1][0][k]
                                *std::pow(L1,n)*std::pow(L2,s-n)
                                /factorial(n)/factorial(s-n);
						}

						for (uint j=32; j<=30+_nf; j++)
						{
							for (uint n=1; n<=s; n++)
                            {
								_B2[j][1][n][k] = recrelNLO_1(_B2[j][0][n-1], k, getExpression("P0ns"));
                                arr.get()[j][k] += _B2[j][1][n][k]*std::pow(L1,n)*std::pow(L2,s-n)/factorial(n)/factorial(s-n);
                            }

                            uint n = 0;
                            double res = recrelNLO_2(_B2[j][0][0], k, 
                                getExpression("P0ns"), 
                                getExpression("P1nsp"));
							_B2[j][1][0][k] = -_B2[j][1][1][k] + res;
                            arr.get()[j][k] += _B2[j][1][0][k]
                                *std::pow(L1,n)*std::pow(L2,s-n)
                                /factorial(n)/factorial(s-n);
						}
					}

                    for (uint j=13; j<=30+_nf; ++j)
                    {
                        for (uint n=0; n<=s; ++n)
                            _B2[j][0][n] = _B2[j][1][n];

                        if (j == (12+_nf))
                            j = 31;
                    }
				}
            } break;
            case 2:
            {
                for (uint j=26; j<=24+_nf; ++j)
                    arr.get()[j] = _C2[j][0][0][0];
                for (uint j=32; j<=30+_nf; ++j)
                    arr.get()[j] = _C2[j][0][0][0];
                arr.get()[25] = _C2[25][0][0][0];

                for (uint s=1; s<_iterations; s++)
				{
					std::println("[DGLAP] NNLO Non-Singlet Iteration {}", s);

					for (uint k=0; k<_grid.size()-1; k++)
					{
						for (uint j=26; j<=24+_nf; j++)
						{
							// recrel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint n=1; n<=t; n++)
                                {
                                    double recrel = recrelNNLO_1(_C2[j][0][t-1][n-1], k, getExpression("P0ns"));
									_C2[j][1][t][n][k] = recrel;

                                    double orig = _C2[j][1][t][n][k];
									double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
									double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
									double res = orig*powers/factorials;

                                    arr.get()[j][k] += res;
                                }
							}

							// recrel #2:
							{
								double fac1 = -0.5*_C2[j][1][s][1][k];
								double fac2 = recrelNNLO_2(_C2[j][0][s-1][0], k, 
                                    getExpression("P0ns"), getExpression("P1nsm"), getExpression("P2nsm"));
								_C2[j][1][s][0][k] = fac1 + fac2;

                                uint n = 0;
                                uint t = s;
                                double orig = _C2[j][1][s][0][k];
                                double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                                double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
                                double res = orig*powers/factorials;

                                arr.get()[j][k] += res;
							}

							// these must be regular ints;
							// unsigned ints, when they are 0 and get --,
							// underflow back to positive 4b,
							// remaining positive and the loop continues (very bad!)

							// recrel #3:
							for (int t=s-1; t>=0; t--)
							{
								double fac1 = -2.0*_alpha_s.beta1()*(_C2[j][1][t+1][0][k] + _C2[j][1][t+1][1][k]);
								double fac2 = recrelNNLO_3(_C2[j][0][t][0], k, 
                                    getExpression("P0ns"), getExpression("P1nsm"));
								_C2[j][1][t][0][k] = fac1 + fac2;

                                uint n = 0;
                                double orig = _C2[j][1][t][0][k];
                                double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                                double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
                                double res = orig*powers/factorials;

                                arr.get()[j][k] += res;
							}
						}

						for (uint j=32; j<=30+_nf; j++)
						{
							// recrel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint n=1; n<=t; n++)
                                {
									double recrel = recrelNNLO_1(_C2[j][0][t-1][n-1], k, getExpression("P0ns"));
									_C2[j][1][t][n][k] = recrel;

                                    double orig = _C2[j][1][t][n][k];
									double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
									double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
									double res = orig*powers/factorials;

                                    arr.get()[j][k] += res;
                                }
							}

							// recrel #2:
							{
								double fac1 = -0.5*_C2[j][1][s][1][k];
								double fac2 = recrelNNLO_2(_C2[j][0][s-1][0], k, 
                                    getExpression("P0ns"), getExpression("P1nsp"), getExpression("P2nsp"));
								_C2[j][1][s][0][k] = fac1 + fac2;

                                uint n = 0;
                                uint t = s;
                                double orig = _C2[j][1][s][0][k];
                                double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                                double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
                                double res = orig*powers/factorials;

                                arr.get()[j][k] += res;
							}

							// these must be regular ints;
							// unsigned ints, when they are 0 and get --,
							// underflow back to positive 4b,
							// remaining positive and the loop continues (very bad!)

							// recrel #3:
							for (int t=s-1; t>=0; t--)
							{
								double fac1 = -2.0*_alpha_s.beta1()*(_C2[j][1][t+1][0][k] + _C2[j][1][t+1][1][k]);
								double fac2 = recrelNNLO_3(_C2[j][0][t][0], k, 
                                    getExpression("P0ns"), getExpression("P1nsp"));
								_C2[j][1][t][0][k] = fac1 + fac2;

                                uint n = 0;
                                double orig = _C2[j][1][t][0][k];
                                double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                                double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
                                double res = orig*powers/factorials;

                                arr.get()[j][k] += res;
							}
						}

					    {
							// recrel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint n=1; n<=t; n++)
                                {
									double recrel = recrelNNLO_1(_C2[25][0][t-1][n-1], k, getExpression("P0ns"));
									_C2[25][1][t][n][k] = recrel;

                                    double orig = _C2[25][1][t][n][k];
									double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
									double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
									double res = orig*powers/factorials;

                                    arr.get()[25][k] += res;
                                }
							}

							// recrel #2:
							{
								double fac1 = -0.5*_C2[25][1][s][1][k];
								double fac2 = recrelNNLO_2(_C2[25][0][s-1][0], k, 
                                    getExpression("P0ns"), getExpression("P1nsm"), getExpression("P2nsv"));
								_C2[25][1][s][0][k] = fac1 + fac2;

                                uint n = 0;
                                uint t = s;
                                double orig = _C2[25][1][s][0][k];
                                double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                                double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
                                double res = orig*powers/factorials;

                                arr.get()[25][k] += res;
							}

							// these must be regular ints;
							// unsigned ints, when they are 0 and get --,
							// underflow back to positive 4b,
							// remaining positive and the loop continues (very bad!)

							// recrel #3:
							for (int t=s-1; t>=0; t--)
							{
								double fac1 = -2.0*_alpha_s.beta1()*(_C2[25][1][t+1][0][k] + _C2[25][1][t+1][1][k]);
								double fac2 = recrelNNLO_3(_C2[25][0][t][0], k, 
                                    getExpression("P0ns"), getExpression("P1nsm"));
								_C2[25][1][t][0][k] = fac1 + fac2;

                                uint n = 0;
                                double orig = _C2[25][1][t][0][k];
                                double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                                double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
                                double res = orig*powers/factorials;

                                arr.get()[25][k] += res;
							}
						}
					}

                    for (uint j=26; j<=24+_nf; j++)
                    {
                        for (uint t=0; t<=s; ++t)
                        {
                            for (uint n=0; n<=t; ++n)
                                _C2[j][0][t][n] = _C2[j][1][t][n];
                        }
                    }
                    for (uint j=32; j<=30+_nf; j++)
                    {
                        for (uint t=0; t<=s; ++t)
                        {
                            for (uint n=0; n<=t; ++n)
                                _C2[j][0][t][n] = _C2[j][1][t][n];
                        }
                    }
                    for (uint t=0; t<=s; ++t)
                    {
                        for (uint n=0; n<=t; ++n)
                            _C2[25][0][t][n] = _C2[25][1][t][n];
                    }
                }
            } break;
            case 3:
            {
                for (uint j=26; j<=24+_nf; ++j)
                    arr.get()[j] = _D2[j][0][0][0][0];
                for (uint j=32; j<=30+_nf; ++j)
                    arr.get()[j] = _D2[j][0][0][0][0];
                arr.get()[25] = _D2[25][0][0][0][0];

                // some shorthand
				double r1 = _r1[_nf];
				double b = _b[_nf];
				double c = _c[_nf];
				double gamma = (r1*r1 + r1*b + c)*_alpha_s.beta3();

                for (uint s=1; s<_iterations; s++)
				{
                    std::println("[DGLAP] N3LO Non-Singlet Iteration {}", s);
					
					for (uint k=0; k<_grid.size()-1; k++)
					{
						// minus distributions
						for (uint j=26; j<=24+_nf; j++)
						{
							// recrel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint m=1; m<=t; m++)
								{
									for (uint n=1; n<=m; n++)
                                    {
										_D2[j][1][t][m][n][k] = recrelN3LO_1(_D2[j][0][t-1][m-1][n-1], k, getExpression("P0ns"));

                                        double orig = _D2[j][1][t][m][n][k];
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

							// recrel #2:
							{
								double fac1 = (
									0.5*(16.0*PI_2*_alpha_s.beta1() + 4*PI*r1*_alpha_s.beta2() - (c + b*r1)*_alpha_s.beta3())
								) * _D2[j][1][s][s][1][k];
								double fac2 = recrelN3LO_2(_D2[j][0][s-1][s-1][0], k,
                                    getExpression("P0ns"), getExpression("P1nsm"), getExpression("P2nsm"), getExpression("P3nsm"));
								_D2[j][1][s][s][0][k] = (fac1 + fac2)/gamma;

                                uint t = s;
                                uint m = s;
                                uint n = 0;
                                double orig = _D2[j][1][s][s][0][k];
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
							for (int m=s-1; m>=0; m--)
							{
								double fac1 = -(
									16.0*PI_2*_alpha_s.beta1() + 4.0*PI*r1*_alpha_s.beta2() + r1*r1*_alpha_s.beta3()
								) * _D2[j][1][s][m+1][1][k];
								double fac2 = recrelN3LO_3(_D2[j][0][s-1][m][0], k,
                                    getExpression("P0ns"), getExpression("P1nsm"), getExpression("P2nsm"), getExpression("P3nsm"));
								_D2[j][1][s][m][0][k] = (fac1 + fac2)/gamma;

                                uint t = s;
                                uint n = 0;
                                double orig = _D2[j][1][s][m][0][k];
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

							// recrel #4:
							for (int t=s-1; t>=0; t--)
							{
								for (int m=t; m>=0; m--)
								{
									double fac1a = -2.0*b*gamma;
									double fac1b = 32*PI_2*(b+r1)*_alpha_s.beta1() - 8*PI*c*_alpha_s.beta2() - 2*c*r1*_alpha_s.beta3();
									double fac1 = fac1a*_D2[j][1][t+1][m+1][0][k] + fac1b*_D2[j][1][t+1][m+1][1][k];
									double fac2 = recrelN3LO_4(_D2[j][0][t][m][0], k,
                                        getExpression("P0ns"), getExpression("P1nsm"), getExpression("P2nsm"), getExpression("P3nsm"));
									_D2[j][1][t][m][0][k] = (fac1 + fac2)/gamma;

                                    uint n = 0;
                                    double orig = _D2[j][1][t][m][0][k];
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

						// plus distributions
						for (uint j=32; j<=30+_nf; j++)
						{
						    // recrel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint m=1; m<=t; m++)
								{
									for (uint n=1; n<=m; n++)
                                    {
										_D2[j][1][t][m][n][k] = recrelN3LO_1(_D2[j][0][t-1][m-1][n-1], k, getExpression("P0ns"));

                                        double orig = _D2[j][1][t][m][n][k];
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

							// recrel #2:
							{
								double fac1 = (
									0.5*(16.0*PI_2*_alpha_s.beta1() + 4*PI*r1*_alpha_s.beta2() - (c + b*r1)*_alpha_s.beta3())
								) * _D2[j][1][s][s][1][k];
								double fac2 = recrelN3LO_2(_D2[j][0][s-1][s-1][0], k,
                                    getExpression("P0ns"), getExpression("P1nsp"), getExpression("P2nsp"), getExpression("P3nsp"));
								_D2[j][1][s][s][0][k] = (fac1 + fac2)/gamma;

                                uint t = s;
                                uint m = s;
                                uint n = 0;
                                double orig = _D2[j][1][s][s][0][k];
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
							for (int m=s-1; m>=0; m--)
							{
								double fac1 = -(
									16.0*PI_2*_alpha_s.beta1() + 4.0*PI*r1*_alpha_s.beta2() + r1*r1*_alpha_s.beta3()
								) * _D2[j][1][s][m+1][1][k];
								double fac2 = recrelN3LO_3(_D2[j][0][s-1][m][0], k,
                                    getExpression("P0ns"), getExpression("P1nsp"), getExpression("P2nsp"), getExpression("P3nsp"));
								_D2[j][1][s][m][0][k] = (fac1 + fac2)/gamma;

                                uint t = s;
                                uint n = 0;
                                double orig = _D2[j][1][s][m][0][k];
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

							// recrel #4:
							for (int t=s-1; t>=0; t--)
							{
								for (int m=t; m>=0; m--)
								{
									double fac1a = -2.0*b*gamma;
									double fac1b = 32*PI_2*(b+r1)*_alpha_s.beta1() - 8*PI*c*_alpha_s.beta2() - 2*c*r1*_alpha_s.beta3();
									double fac1 = fac1a*_D2[j][1][t+1][m+1][0][k] + fac1b*_D2[j][1][t+1][m+1][1][k];
									double fac2 = recrelN3LO_4(_D2[j][0][t][m][0], k,
                                        getExpression("P0ns"), getExpression("P1nsp"), getExpression("P2nsp"), getExpression("P3nsp"));
									_D2[j][1][t][m][0][k] = (fac1 + fac2)/gamma;

                                    uint n = 0;
                                    double orig = _D2[j][1][t][m][0][k];
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

						// valence distribution
						{
						    // recrel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint m=1; m<=t; m++)
								{
									for (uint n=1; n<=m; n++)
                                    {
										_D2[25][1][t][m][n][k] = recrelN3LO_1(_D2[25][0][t-1][m-1][n-1], k, getExpression("P0ns"));

                                        double orig = _D2[25][1][t][m][n][k];
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
										
										arr.get()[25][k] += res;
                                    }
								}
							}

							// recrel #2:
							{
								double fac1 = (
									0.5*(16.0*PI_2*_alpha_s.beta1() + 4*PI*r1*_alpha_s.beta2() - (c + b*r1)*_alpha_s.beta3())
								) * _D2[25][1][s][s][1][k];
								double fac2 = recrelN3LO_2(_D2[25][0][s-1][s-1][0], k,
                                    getExpression("P0ns"), getExpression("P1nsm"), getExpression("P2nsv"), getExpression("P3nsv"));
								_D2[25][1][s][s][0][k] = (fac1 + fac2)/gamma;

                                uint t = s;
                                uint m = s;
                                uint n = 0;
                                double orig = _D2[25][1][s][s][0][k];
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
                                
                                arr.get()[25][k] += res;
							}

							// recrel #3:
							for (int m=s-1; m>=0; m--)
							{
								double fac1 = -(
									16.0*PI_2*_alpha_s.beta1() + 4.0*PI*r1*_alpha_s.beta2() + r1*r1*_alpha_s.beta3()
								) * _D2[25][1][s][m+1][1][k];
								double fac2 = recrelN3LO_3(_D2[25][0][s-1][m][0], k,
                                    getExpression("P0ns"), getExpression("P1nsm"), getExpression("P2nsv"), getExpression("P3nsv"));
								_D2[25][1][s][m][0][k] = (fac1 + fac2)/gamma;

                                uint t = s;
                                uint n = 0;
                                double orig = _D2[25][1][s][m][0][k];
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
                                
                                arr.get()[25][k] += res;
							}

							// recrel #4:
							for (int t=s-1; t>=0; t--)
							{
								for (int m=t; m>=0; m--)
								{
									double fac1a = -2.0*b*gamma;
									double fac1b = 32*PI_2*(b+r1)*_alpha_s.beta1() - 8*PI*c*_alpha_s.beta2() - 2*c*r1*_alpha_s.beta3();
									double fac1 = fac1a*_D2[25][1][t+1][m+1][0][k] + fac1b*_D2[25][1][t+1][m+1][1][k];
									double fac2 = recrelN3LO_4(_D2[25][0][t][m][0], k,
                                        getExpression("P0ns"), getExpression("P1nsm"), getExpression("P2nsv"), getExpression("P3nsv"));
									_D2[25][1][t][m][0][k] = (fac1 + fac2)/gamma;

                                    uint n = 0;
                                    double orig = _D2[25][1][t][m][0][k];
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
                                    
                                    arr.get()[25][k] += res;
								}
							}
						}
					}

                    for (uint j=26; j<=24+_nf; j++)
                    {
                        for (uint t=0; t<=s; ++t)
                        {
                            for (uint m=0; m<=s; ++m)
                            {
                                for (uint n=0; n<=s; ++n)
                                    _D2[j][0][t][m][n] = _D2[j][1][t][m][n];
                            }
                        }
                    }
                    for (uint j=32; j<=30+_nf; j++)
                    {
                        for (uint t=0; t<=s; ++t)
                        {
                            for (uint m=0; m<=s; ++m)
                            {
                                for (uint n=0; n<=s; ++n)
                                    _D2[j][0][t][m][n] = _D2[j][1][t][m][n];
                            }
                        }
                    }
                    for (uint t=0; t<=s; ++t)
                    {
                        for (uint m=0; m<=s; ++m)
                        {
                            for (uint n=0; n<=s; ++n)
                                _D2[25][0][t][m][n] = _D2[25][1][t][m][n];
                        }
                    }
                }
            } break;
        }
    }
}
