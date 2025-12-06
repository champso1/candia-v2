#include "Candia-v2/Candia.hpp"
#include "Candia-v2/FuncArrGrid.hpp"
#include "Candia-v2/Math.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"

#include <functional>
#include <print>
#include <thread>

namespace Candia2
{
    void DGLAPSolver::loadAllExpressions()
    {
        createExpressionGrid<P0ns>("P0ns", _grid);
        createExpressionGrid<P0qq>("P0qq", _grid);
        createExpressionGrid<P0qg>("P0qg", _grid);
        createExpressionGrid<P0gq>("P0gq", _grid);
        createExpressionGrid<P0gg>("P0gg", _grid);
    
        if (_order >= 1)
        {
            createExpressionGrid<P1nsm>("P1nsm", _grid);
            createExpressionGrid<P1nsp>("P1nsp", _grid);
            createExpressionGrid<P1qq>("P1qq", _grid);
            createExpressionGrid<P1qg>("P1qg", _grid);
            createExpressionGrid<P1gq>("P1gq", _grid);
            createExpressionGrid<P1gg>("P1gg", _grid);
        }
        if (_order >= 2)
        {
            createExpressionGrid<P2nsm>("P2nsm", _grid);
            createExpressionGrid<P2nsp>("P2nsp", _grid);
            createExpressionGrid<P2nsv>("P2nsv", _grid);
            createExpressionGrid<P2qq>("P2qq", _grid);
            createExpressionGrid<P2qg>("P2qg", _grid);
            createExpressionGrid<P2gq>("P2gq", _grid);
            createExpressionGrid<P2gg>("P2gg", _grid);

            createExpressionGrid<A2ns>("A2ns", _grid);
            createExpressionGrid<A2gq>("A2gq", _grid);
            createExpressionGrid<A2gg>("A2gg", _grid);
            createExpressionGrid<A2hq>("A2hq", _grid);
            createExpressionGrid<A2hg>("A2hg", _grid);
        }
        if (_order >= 3)
        {
            createExpressionGrid<P3nsm>("P3nsm", _grid);
            createExpressionGrid<P3nsp>("P3nsp", _grid);
            createExpressionGrid<P3nsv>("P3nsv", _grid);
            createExpressionGrid<P3qq>("P3qq", _grid);
            createExpressionGrid<P3qg>("P3qg", _grid);
            createExpressionGrid<P3gq>("P3gq", _grid);
            createExpressionGrid<P3gg>("P3gg", _grid);

            
            createExpressionGrid<OpMatElemN3LO>("A3nsm",  _grid, ome::AqqQNSEven);
            createExpressionGrid<OpMatElemN3LO>("A3nsp",  _grid, ome::AqqQNSOdd);
            createExpressionGrid<OpMatElemN3LO>("A3gq",   _grid, ome::AgqQ);
            createExpressionGrid<OpMatElemN3LO>("A3gg",   _grid, ome::AggQ);
            createExpressionGrid<OpMatElemN3LO>("A3hq",   _grid, ome::AQqPS);
            createExpressionGrid<OpMatElemN3LO>("A3hg",   _grid, ome::AQg);
            createExpressionGrid<OpMatElemN3LO>("A3psqq", _grid, ome::AqqQPS);
            createExpressionGrid<OpMatElemN3LO>("A3sqg",  _grid, ome::AqgQ);
        }
        
    }

    void DGLAPSolver::SetupCoefficients2()
    {
        switch (_order)
		{
			case 0: // LO
			{
				for (uint k=0; k<_grid.Size(); k++)
                {
                    for (uint j=13; j<=18; j++)
                        _A2[j][0][k] = _A2[j-12][0][k]-_A2[j-6][0][k];
                
                    _A2[25][0][k]=0.;
                    for (uint j=13; j<=18; j++)
                        _A2[25][0][k] += _A2[j][0][k];

                    for (uint j=26; j<=30; j++)
                        _A2[j][0][k] = _A2[13][0][k]-_A2[j-12][0][k];

                    for (uint j=19; j<=24; j++)
                        _A2[j][0][k] = _A2[j-18][0][k]+_A2[j-12][0][k];

                    _S2[0][1][0][k] = 0.0;
                    for (uint j=19; j<=24; j++)
                        _S2[0][1][0][k] += _A2[j][0][k];
                
                    for (uint j=32; j<=36; j++)
                        _A2[j][0][k]=_A2[19][0][k]-_A2[j-12][0][k];
                }
			} break;
			case 1: // NLO
			{
				for (uint k=0; k<_grid.Size(); k++)
				{
					for (uint j=13; j<=18; j++)
						_B2[j][0][0][k] = _B2[j-12][0][0][k]-_B2[j-6][0][0][k];
			
					_B2[25][0][0][k]=0.;
					for (uint j=13; j<=18; j++)
						_B2[25][0][0][k] += _B2[j][0][0][k];
					for (uint j=26; j<=30; j++)
						_B2[j][0][0][k] = _B2[13][0][0][k]-_B2[j-12][0][0][k];
					for (uint j=19; j<=24; j++)
						_B2[j][0][0][k] = _B2[j-18][0][0][k]+_B2[j-12][0][0][k];

					_S2[0][1][0][k] = 0.0;
					for (uint j=19; j<=24; j++)
						_S2[0][1][0][k] += _B2[j][0][0][k];
			
					for (uint j=32; j<=36; j++)
						_B2[j][0][0][k] = _B2[19][0][0][k]-_B2[j-12][0][0][k];
				}
			} break;
			case 2: // NNLO
			{
				for (uint k=0; k<_grid.Size(); k++)
				{
					for (uint j=13; j<=18; j++)
						_C2[j][0][0][0][k] = _C2[j-12][0][0][0][k]-_C2[j-6][0][0][0][k];
			
					_C2[25][0][0][0][k]=0.0;
					for (uint j=13; j<=18; j++)
						_C2[25][0][0][0][k] += _C2[j][0][0][0][k];

					for (uint j=26; j<=30; j++)
						_C2[j][0][0][0][k] = _C2[13][0][0][0][k]-_C2[j-12][0][0][0][k];

					for (uint j=19; j<=24; j++)
						_C2[j][0][0][0][k] = _C2[j-18][0][0][0][k]+_C2[j-12][0][0][0][k];

					_S2[0][1][0][k] = 0.0;
					for (uint j=19; j<=24; j++)
						_S2[0][1][0][k] += _C2[j][0][0][0][k];
			
					for (uint j=32; j<=36; j++)
						_C2[j][0][0][0][k] = _C2[19][0][0][0][k]-_C2[j-12][0][0][0][k];
				}
			} break;
			case 3: //N3LO
			{
				for (uint k=0; k<_grid.Size(); k++)
				{
					for (uint j=13; j<=18; j++)
						_D2[j][0][0][0][0][k] = _D2[j-12][0][0][0][0][k]-_D2[j-6][0][0][0][0][k];
			
					_D2[25][0][0][0][0][k]=0.;
					for (uint j=13; j<=18; j++)
						_D2[25][0][0][0][0][k] += _D2[j][0][0][0][0][k];
					for (uint j=26; j<=30; j++)
						_D2[j][0][0][0][0][k] = _D2[13][0][0][0][0][k]-_D2[j-12][0][0][0][0][k];
					for (uint j=19; j<=24; j++)
						_D2[j][0][0][0][0][k] = _D2[j-18][0][0][0][0][k]+_D2[j-12][0][0][0][0][k];

					_S2[0][1][0][k] = 0.0;
					for (uint j=19; j<=24; j++)
						_S2[0][1][0][k] += _D2[j][0][0][0][0][k];
			
					for (uint j=32; j<=36; j++)
						_D2[j][0][0][0][0][k] = _D2[19][0][0][0][0][k]-_D2[j-12][0][0][0][0][k];
				}
			}
		}
    }


    void DGLAPSolver::EvolveSinglet2(
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
                    std::cerr << "[DGLAP2] LO Singlet Iteration " << n << '\n';
                    
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        _S2[0][1][1][k] = RecRelS_1(_S2[0][1][0], k, getSplitFunc("P0qq")) + RecRelS_1(_S2[0][0][0], k, getSplitFunc("P0qg"));
                        _S2[0][0][1][k] = RecRelS_1(_S2[0][1][0], k, getSplitFunc("P0gq")) + RecRelS_1(_S2[0][0][0], k, getSplitFunc("P0gg"));

                        for (uint t=0; t<=_trunc_idx; ++t)
                        {
                            for (uint j=0; j<=1; j++)
                                arr.get()[j*31][k] += _S2[t][j][1][k] * std::pow(_alpha1, t) * std::pow(L1, n)/Factorial(n);
                        }

                    }

                    for (uint j=0; j<=1; ++j)
                        _S2[0][j][0] = _S2[0][j][1];
                }
            } break;
            case 1:
            {
                for (uint n=1; n<_iterations; n++)
                {
                    std::cerr << "[DGLAP2] NLO Singlet Iteration " << n << '\n';
                    
                    // LO piece (non truncated)
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        _S2[0][1][1][k] = RecRelS_1(_S2[0][1][0], k, getSplitFunc("P0qq")) + RecRelS_1(_S2[0][0][0], k, getSplitFunc("P0qg"));
                        _S2[0][0][1][k] = RecRelS_1(_S2[0][1][0], k, getSplitFunc("P0gq")) + RecRelS_1(_S2[0][0][0], k, getSplitFunc("P0gg"));
                    }

                    // new NLO piece non-convolution
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        for (uint j=0; j<=1; j++)
                            _S2[1][j][1][k] = -_S2[0][j][1][k] * _alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()) - _S2[1][j][0][k];
                    }

                    // new NLO piece convolution
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        _S2[1][1][1][k]+=RecRelS_2(_S2[0][1][0],_S2[1][1][0],k,getSplitFunc("P0qq"),getSplitFunc("P1qq"))+RecRelS_2(_S2[0][0][0],_S2[1][0][0],k,getSplitFunc("P0qg"),getSplitFunc("P1qg"));
					    _S2[1][0][1][k]+=RecRelS_2(_S2[0][1][0],_S2[1][1][0],k,getSplitFunc("P0gq"),getSplitFunc("P1gq"))+RecRelS_2(_S2[0][0][0],_S2[1][0][0],k,getSplitFunc("P0gg"),getSplitFunc("P1gg"));
                    }

                    // NLO truncation terms
                    for (uint t=2; t<=_trunc_idx; ++t)
                    {
                        double T = static_cast<double>(t);

                        // non-convolution piece:
                        for (uint k=0; k<_grid.Size()-1; k++)
                        {
                            for (uint j=0; j<=1; j++)
                                _S2[t][j][1][k] =
                                    - (_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()))*_S2[t-1][j][1][k]
                                    - T*_S2[t][j][0][k]
                                    - (T-1.0)*(_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()))*_S2[t-1][j][0][k];
                        }

                        // convolution piece
                        for (uint k=0; k<_grid.Size()-1; k++)
                        {
                            _S2[t][1][1][k]+=RecRelS_2(_S2[t-1][1][0],_S2[t][1][0],k,getSplitFunc("P0qq"),getSplitFunc("P1qq"))+RecRelS_2(_S2[t-1][0][0],_S2[t][0][0],k,getSplitFunc("P0qg"),getSplitFunc("P1qg"));
                            _S2[t][0][1][k]+=RecRelS_2(_S2[t-1][1][0],_S2[t][1][0],k,getSplitFunc("P0gq"),getSplitFunc("P1gq"))+RecRelS_2(_S2[t-1][0][0],_S2[t][0][0],k,getSplitFunc("P0gg"),getSplitFunc("P1gg"));
                        }
                    }

                    for (uint t=0; t<=_trunc_idx; ++t)
                    {
                        for (uint j=0; j<=1; ++j)
                        {
                            for (uint k=0; k<_grid.Size()-1; k++)
                                arr.get()[j*31][k] += _S2[t][j][1][k] * std::pow(_alpha1, t) * std::pow(L1, n)/Factorial(n);        
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
                    std::cerr << "[DGLAP2] NNLO Singlet Iteration " << n << '\n';

                    // LO piece (non truncated)
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        _S2[0][1][1][k] = RecRelS_1(_S2[0][1][0], k, getSplitFunc("P0qq")) + RecRelS_1(_S2[0][0][0], k, getSplitFunc("P0qg"));
                        _S2[0][0][1][k] = RecRelS_1(_S2[0][1][0], k, getSplitFunc("P0gq")) + RecRelS_1(_S2[0][0][0], k, getSplitFunc("P0gg"));
                    }

                    // new NLO piece non-convolution
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        for (uint j=0; j<=1; j++)
                            _S2[1][j][1][k] = -_S2[0][j][1][k] * _alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()) - _S2[1][j][0][k];
                    }

                    // new NLO piece convolution
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        _S2[1][1][1][k]+=RecRelS_2(_S2[0][1][0],_S2[1][1][0],k,getSplitFunc("P0qq"),getSplitFunc("P1qq"))+RecRelS_2(_S2[0][0][0],_S2[1][0][0],k,getSplitFunc("P0qg"),getSplitFunc("P1qg"));
					    _S2[1][0][1][k]+=RecRelS_2(_S2[0][1][0],_S2[1][1][0],k,getSplitFunc("P0gq"),getSplitFunc("P1gq"))+RecRelS_2(_S2[0][0][0],_S2[1][0][0],k,getSplitFunc("P0gg"),getSplitFunc("P1gg"));
                    }

                    // new NNLO piece non-convolution
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        for (uint j=0; j<=1; j++)
                            _S2[2][j][1][k] =
                                - (_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()))*_S2[1][j][1][k]
                                - (_alpha_s.Beta2()/(16.0*PI_2*_alpha_s.Beta0()))*_S2[0][j][1][k]
                                - 2.0*_S2[2][j][0][k]
                                - (_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()))*_S2[1][j][0][k];
                    }

                    // new NNLO piece non-convolution
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        _S2[2][1][1][k] += 
                            RecRelS_3(
                                _S2[0][1][0], _S2[1][1][0], _S2[2][1][0], k, 
                                getSplitFunc("P0qq"), getSplitFunc("P1qq"), getSplitFunc("P2qq")) +
                            RecRelS_3(
                                _S2[0][0][0], _S2[1][0][0], _S2[2][0][0], k,
                                getSplitFunc("P0qg"), getSplitFunc("P1qg"), getSplitFunc("P2qg"));

                        _S2[2][0][1][k] += 
                            RecRelS_3(_S2[0][1][0], _S2[1][1][0], _S2[2][1][0], k,
                                getSplitFunc("P0gq"), getSplitFunc("P1gq"), getSplitFunc("P2gq")) +
                            RecRelS_3(_S2[0][0][0], _S2[1][0][0], _S2[2][0][0], k,
                                getSplitFunc("P0gg"), getSplitFunc("P1gg"), getSplitFunc("P2gg"));
                    }

                    for (uint t=3; t<=_trunc_idx; ++t)
                    {
                        double T = static_cast<double>(t);

                        // non-convolution piece:
                        for (uint k=0; k<_grid.Size()-1; k++)
                        {
                            for (uint j=0; j<=1; j++)
                                _S2[t][j][1][k] =
                                    - (_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()))*_S2[t-1][j][1][k]
                                    - (_alpha_s.Beta2()/(16.0*PI_2*_alpha_s.Beta0()))*_S2[t-2][j][1][k]
                                    - T*_S2[t][j][0][k]
                                    - (T-1.0)*(_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()))*_S2[t-1][j][0][k]
                                    - (T-2.0)*(_alpha_s.Beta2()/(16.0*PI_2*_alpha_s.Beta0()))*_S2[t-2][j][0][k];
                        }

                        // convolution piece
                        for (uint k=0; k<_grid.Size()-1; k++)
                        {
                            _S2[t][1][1][k] += 
                                RecRelS_3(_S2[t-2][1][0], _S2[t-1][1][0], _S2[t][1][0], k,
                                    getSplitFunc("P0qq"), getSplitFunc("P1qq"), getSplitFunc("P2qq")) +
                                RecRelS_3(_S2[t-2][0][0], _S2[t-1][0][0], _S2[t][0][0], k,
                                    getSplitFunc("P0qg"), getSplitFunc("P1qg"), getSplitFunc("P2qg"));

                            _S2[t][0][1][k] += 
                                RecRelS_3(_S2[t-2][1][0], _S2[t-1][1][0], _S2[t][1][0], k,
                                    getSplitFunc("P0gq"), getSplitFunc("P1gq"), getSplitFunc("P2gq"))+
                                RecRelS_3(_S2[t-2][0][0], _S2[t-1][0][0], _S2[t][0][0], k,
                                    getSplitFunc("P0gg"), getSplitFunc("P1gg"), getSplitFunc("P2gg"));
                        }
                    }

                    for (uint t=0; t<=_trunc_idx; ++t)
                    {
                        for (uint j=0; j<=1; ++j)
                        {
                            for (uint k=0; k<_grid.Size()-1; k++)
                                arr.get()[j*31][k] += _S2[t][j][1][k] * std::pow(_alpha1, t) * std::pow(L1, n)/Factorial(n);        
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
                    std::cerr << "[DGLAP2] N3LO Singlet Iteration " << n << '\n';

                    // LO piece
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        _S2[0][1][1][k] = RecRelS_1(_S2[0][1][0], k, getSplitFunc("P0qq")) + RecRelS_1(_S2[0][0][0], k, getSplitFunc("P0qg"));
                        _S2[0][0][1][k] = RecRelS_1(_S2[0][1][0], k, getSplitFunc("P0gq")) + RecRelS_1(_S2[0][0][0], k, getSplitFunc("P0gg"));
                    }

                    // new NLO piece non-convolution
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        for (uint j=0; j<=1; j++)
                            _S2[1][j][1][k] = -_S2[0][j][1][k] * _alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()) - _S2[1][j][0][k];
                    }

                    // new NLO piece convolution
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        _S2[1][1][1][k]+=RecRelS_2(_S2[0][1][0],_S2[1][1][0],k,getSplitFunc("P0qq"),getSplitFunc("P1qq"))+RecRelS_2(_S2[0][0][0],_S2[1][0][0],k,getSplitFunc("P0qg"),getSplitFunc("P1qg"));
					    _S2[1][0][1][k]+=RecRelS_2(_S2[0][1][0],_S2[1][1][0],k,getSplitFunc("P0gq"),getSplitFunc("P1gq"))+RecRelS_2(_S2[0][0][0],_S2[1][0][0],k,getSplitFunc("P0gg"),getSplitFunc("P1gg"));
                    }

                    // new NNLO piece non-convolution
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        for (uint j=0; j<=1; j++)
                            _S2[2][j][1][k] =
                                - (_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()))*_S2[1][j][1][k]
                                - (_alpha_s.Beta2()/(16.0*PI_2*_alpha_s.Beta0()))*_S2[0][j][1][k]
                                - 2.0*_S2[2][j][0][k]
                                - (_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()))*_S2[1][j][0][k];
                    }

                    // new NNLO piece non-convolution
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        _S2[2][1][1][k] += 
                            RecRelS_3(
                                _S2[0][1][0], _S2[1][1][0], _S2[2][1][0], k, 
                                getSplitFunc("P0qq"), getSplitFunc("P1qq"), getSplitFunc("P2qq")) +
                            RecRelS_3(
                                _S2[0][0][0], _S2[1][0][0], _S2[2][0][0], k,
                                getSplitFunc("P0qg"), getSplitFunc("P1qg"), getSplitFunc("P2qg"));

                        _S2[2][0][1][k] += 
                            RecRelS_3(_S2[0][1][0], _S2[1][1][0], _S2[2][1][0], k,
                                getSplitFunc("P0gq"), getSplitFunc("P1gq"), getSplitFunc("P2gq")) +
                            RecRelS_3(_S2[0][0][0], _S2[1][0][0], _S2[2][0][0], k,
                                getSplitFunc("P0gg"), getSplitFunc("P1gg"), getSplitFunc("P2gg"));
                    }

                    // new N3LO piece non-convolution
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        for (uint j=0; j<=1; j++)
                            _S2[3][j][1][k] =
                                - (_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()))*_S2[2][j][1][k]
                                - (_alpha_s.Beta2()/(16.0*PI_2*_alpha_s.Beta0()))*_S2[1][j][1][k]
                                - (_alpha_s.Beta3()/(64.0*PI_3*_alpha_s.Beta0()))*_S2[0][j][1][k]
                                - 3.0*_S2[3][j][0][k]
                                - 2.0*(_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()))*_S2[2][j][0][k]
                                - (_alpha_s.Beta2()/(16.0*PI_2*_alpha_s.Beta0()))*_S2[1][j][0][k];
                    }

                    // new N3LO piece convolution
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        _S2[3][1][1][k] += 
                            RecRelS_4(_S2[0][1][0], _S2[1][1][0], _S2[2][1][0], _S2[3][1][0], k,
                                getSplitFunc("P0qq"), getSplitFunc("P1qq"), getSplitFunc("P2qq"), getSplitFunc("P3qq")) +
                            RecRelS_4(_S2[0][0][0], _S2[1][0][0], _S2[2][0][0], _S2[3][0][0], k,
                                getSplitFunc("P0qg"), getSplitFunc("P1qg"), getSplitFunc("P2qg"), getSplitFunc("P3qg"));

                        _S2[3][0][1][k] += 
                            RecRelS_4(_S2[0][1][0], _S2[1][1][0], _S2[2][1][0], _S2[3][1][0], k,
                                getSplitFunc("P0gq"), getSplitFunc("P1gq"), getSplitFunc("P2gq"), getSplitFunc("P3gq")) +
                            RecRelS_4(_S2[0][0][0], _S2[1][0][0], _S2[2][0][0], _S2[3][0][0], k,
                                getSplitFunc("P0gg"), getSplitFunc("P1gg"), getSplitFunc("P2gg"), getSplitFunc("P3gg"));
                    }

                    for (uint t=4; t<=_trunc_idx; ++t)
                    {
                        double T = static_cast<double>(t);

                        // non-convolution piece:
                        for (uint k=0; k<_grid.Size()-1; k++)
                        {
                            for (uint j=0; j<=1; j++)
                                _S2[t][j][1][k] =
                                    - (_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()))*_S2[t-1][j][1][k]
                                    - (_alpha_s.Beta2()/(16.0*PI_2*_alpha_s.Beta0()))*_S2[t-2][j][1][k]
                                    - (_alpha_s.Beta3()/(64.0*PI_3*_alpha_s.Beta0()))*_S2[t-3][j][1][k]
                                    - T*_S2[t][j][0][k]
                                    - (T-1.0)*(_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()))*_S2[t-1][j][0][k]
                                    - (T-2.0)*(_alpha_s.Beta2()/(16.0*PI_2*_alpha_s.Beta0()))*_S2[t-2][j][0][k]
                                    - (T-3.0)*(_alpha_s.Beta3()/(64.0*PI_3*_alpha_s.Beta0()))*_S2[t-3][j][0][k];
                        }

                        // convolution piece
                        for (uint k=0; k<_grid.Size()-1; k++)
                        {
                            _S2[t][1][1][k] += 
                                RecRelS_4(_S2[t-3][1][0], _S2[t-2][1][0], _S2[t-1][1][0], _S2[t][1][0],k,
                                    getSplitFunc("P0qq"), getSplitFunc("P1qq"), getSplitFunc("P2qq"), getSplitFunc("P3qq")) +
                                RecRelS_4(_S2[t-3][0][0], _S2[t-2][0][0], _S2[t-1][0][0], _S2[t][0][0],k,
                                    getSplitFunc("P0qg"), getSplitFunc("P1qg"), getSplitFunc("P2qg"), getSplitFunc("P3qg"));
                            _S2[t][0][1][k] += 
                                RecRelS_4(_S2[t-3][1][0], _S2[t-2][1][0], _S2[t-1][1][0], _S2[t][1][0],k,
                                    getSplitFunc("P0gq"), getSplitFunc("P1gq"), getSplitFunc("P2gq"), getSplitFunc("P3gq")) +
                                RecRelS_4(_S2[t-3][0][0], _S2[t-2][0][0], _S2[t-1][0][0], _S2[t][0][0],k,
                                    getSplitFunc("P0gg"), getSplitFunc("P1gg"), getSplitFunc("P2gg"), getSplitFunc("P3gg"));
                        }
                    }

                    for (uint t=0; t<=_trunc_idx; ++t)
                    {
                        for (uint j=0; j<=1; ++j)
                        {
                            for (uint k=0; k<_grid.Size()-1; k++)
                                arr.get()[j*31][k] += _S2[t][j][1][k] * std::pow(_alpha1, t) * std::pow(L1, n)/Factorial(n);        
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



    void DGLAPSolver::EvolveNonSinglet2(
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
                    std::cerr << "[DGLAP2] LO Non-Singlet Iteration " << n << '\n';
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        for (uint j=13; j<=30+_nf; j++)
                        {
                            _A2[j][1][k] = RecRelLO(_A2[j][0], k, getSplitFunc("P0ns"));
                            arr.get()[j][k] += _A2[j][1][k]*std::pow(L1, n)/Factorial(n);

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
					std::cerr << "[DGLAP2] NLO Non-Singlet Iteration " << s << '\n';

					for (uint k=0; k<_grid.Size()-1;k++)
					{
						for (uint j=13; j<=12+_nf; j++)
						{
							for (uint n=1; n<=s; n++)
                            {
								_B2[j][1][n][k] = RecRelNLO_1(_B2[j][0][n-1], k, getSplitFunc("P0ns"));
                                arr.get()[j][k] += _B2[j][1][n][k]*std::pow(L1,n)*std::pow(L2,s-n)/Factorial(n)/Factorial(s-n);
                            }

                            uint n = 0;
                            double res = RecRelNLO_2(_B2[j][0][0], k, 
                                getSplitFunc("P0ns"), 
                                getSplitFunc("P1nsm"));
							_B2[j][1][0][k] = -_B2[j][1][1][k] + res;
                            arr.get()[j][k] += _B2[j][1][0][k]
                                *std::pow(L1,n)*std::pow(L2,s-n)
                                /Factorial(n)/Factorial(s-n);
						}

						for (uint j=32; j<=30+_nf; j++)
						{
							for (uint n=1; n<=s; n++)
                            {
								_B2[j][1][n][k] = RecRelNLO_1(_B2[j][0][n-1], k, getSplitFunc("P0ns"));
                                arr.get()[j][k] += _B2[j][1][n][k]*std::pow(L1,n)*std::pow(L2,s-n)/Factorial(n)/Factorial(s-n);
                            }

                            uint n = 0;
                            double res = RecRelNLO_2(_B2[j][0][0], k, 
                                getSplitFunc("P0ns"), 
                                getSplitFunc("P1nsp"));
							_B2[j][1][0][k] = -_B2[j][1][1][k] + res;
                            arr.get()[j][k] += _B2[j][1][0][k]
                                *std::pow(L1,n)*std::pow(L2,s-n)
                                /Factorial(n)/Factorial(s-n);
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
					std::cerr << "[DGLAP2] NNLO Non-Singlet Iteration " << s << '\n';

					for (uint k=0; k<_grid.Size()-1; k++)
					{
						for (uint j=26; j<=24+_nf; j++)
						{
							// RecRel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint n=1; n<=t; n++)
                                {
                                    double recrel = RecRelNNLO_1(_C2[j][0][t-1][n-1], k, getSplitFunc("P0ns"));
									_C2[j][1][t][n][k] = recrel;

                                    double orig = _C2[j][1][t][n][k];
									double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
									double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
									double res = orig*powers/factorials;

                                    arr.get()[j][k] += res;
                                }
							}

							// RecRel #2:
							{
								double fac1 = -0.5*_C2[j][1][s][1][k];
								double fac2 = RecRelNNLO_2(_C2[j][0][s-1][0], k, 
                                    getSplitFunc("P0ns"), getSplitFunc("P1nsm"), getSplitFunc("P2nsm"));
								_C2[j][1][s][0][k] = fac1 + fac2;

                                uint n = 0;
                                uint t = s;
                                double orig = _C2[j][1][s][0][k];
                                double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                                double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
                                double res = orig*powers/factorials;

                                arr.get()[j][k] += res;
							}

							// these must be regular ints;
							// unsigned ints, when they are 0 and get --,
							// underflow back to positive 4b,
							// remaining positive and the loop continues (very bad!)

							// RecRel #3:
							for (int t=s-1; t>=0; t--)
							{
								double fac1 = -2.0*_alpha_s.Beta1()*(_C2[j][1][t+1][0][k] + _C2[j][1][t+1][1][k]);
								double fac2 = RecRelNNLO_3(_C2[j][0][t][0], k, 
                                    getSplitFunc("P0ns"), getSplitFunc("P1nsm"));
								_C2[j][1][t][0][k] = fac1 + fac2;

                                uint n = 0;
                                double orig = _C2[j][1][t][0][k];
                                double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                                double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
                                double res = orig*powers/factorials;

                                arr.get()[j][k] += res;
							}
						}

						for (uint j=32; j<=30+_nf; j++)
						{
							// RecRel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint n=1; n<=t; n++)
                                {
									double recrel = RecRelNNLO_1(_C2[j][0][t-1][n-1], k, getSplitFunc("P0ns"));
									_C2[j][1][t][n][k] = recrel;

                                    double orig = _C2[j][1][t][n][k];
									double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
									double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
									double res = orig*powers/factorials;

                                    arr.get()[j][k] += res;
                                }
							}

							// RecRel #2:
							{
								double fac1 = -0.5*_C2[j][1][s][1][k];
								double fac2 = RecRelNNLO_2(_C2[j][0][s-1][0], k, 
                                    getSplitFunc("P0ns"), getSplitFunc("P1nsp"), getSplitFunc("P2nsp"));
								_C2[j][1][s][0][k] = fac1 + fac2;

                                uint n = 0;
                                uint t = s;
                                double orig = _C2[j][1][s][0][k];
                                double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                                double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
                                double res = orig*powers/factorials;

                                arr.get()[j][k] += res;
							}

							// these must be regular ints;
							// unsigned ints, when they are 0 and get --,
							// underflow back to positive 4b,
							// remaining positive and the loop continues (very bad!)

							// RecRel #3:
							for (int t=s-1; t>=0; t--)
							{
								double fac1 = -2.0*_alpha_s.Beta1()*(_C2[j][1][t+1][0][k] + _C2[j][1][t+1][1][k]);
								double fac2 = RecRelNNLO_3(_C2[j][0][t][0], k, 
                                    getSplitFunc("P0ns"), getSplitFunc("P1nsp"));
								_C2[j][1][t][0][k] = fac1 + fac2;

                                uint n = 0;
                                double orig = _C2[j][1][t][0][k];
                                double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                                double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
                                double res = orig*powers/factorials;

                                arr.get()[j][k] += res;
							}
						}

					    {
							// RecRel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint n=1; n<=t; n++)
                                {
									double recrel = RecRelNNLO_1(_C2[25][0][t-1][n-1], k, getSplitFunc("P0ns"));
									_C2[25][1][t][n][k] = recrel;

                                    double orig = _C2[25][1][t][n][k];
									double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
									double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
									double res = orig*powers/factorials;

                                    arr.get()[25][k] += res;
                                }
							}

							// RecRel #2:
							{
								double fac1 = -0.5*_C2[25][1][s][1][k];
								double fac2 = RecRelNNLO_2(_C2[25][0][s-1][0], k, 
                                    getSplitFunc("P0ns"), getSplitFunc("P1nsm"), getSplitFunc("P2nsv"));
								_C2[25][1][s][0][k] = fac1 + fac2;

                                uint n = 0;
                                uint t = s;
                                double orig = _C2[25][1][s][0][k];
                                double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                                double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
                                double res = orig*powers/factorials;

                                arr.get()[25][k] += res;
							}

							// these must be regular ints;
							// unsigned ints, when they are 0 and get --,
							// underflow back to positive 4b,
							// remaining positive and the loop continues (very bad!)

							// RecRel #3:
							for (int t=s-1; t>=0; t--)
							{
								double fac1 = -2.0*_alpha_s.Beta1()*(_C2[25][1][t+1][0][k] + _C2[25][1][t+1][1][k]);
								double fac2 = RecRelNNLO_3(_C2[25][0][t][0], k, 
                                    getSplitFunc("P0ns"), getSplitFunc("P1nsm"));
								_C2[25][1][t][0][k] = fac1 + fac2;

                                uint n = 0;
                                double orig = _C2[25][1][t][0][k];
                                double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                                double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
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
				double gamma = (r1*r1 + r1*b + c)*_alpha_s.Beta3();

                for (uint s=1; s<_iterations; s++)
				{
                    std::cerr << "[DGLAP2] N3LO Non-Singlet Iteration " << s << '\n';
					
					for (uint k=0; k<_grid.Size()-1; k++)
					{
						// minus distributions
						for (uint j=26; j<=24+_nf; j++)
						{
							// RecRel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint m=1; m<=t; m++)
								{
									for (uint n=1; n<=m; n++)
                                    {
										_D2[j][1][t][m][n][k] = RecRelN3LO_1(_D2[j][0][t-1][m-1][n-1], k, getSplitFunc("P0ns"));

                                        double orig = _D2[j][1][t][m][n][k];
										double powers =
											std::pow(L1,n)
											*std::pow(L2,(m-n))
											*std::pow(L3,(t-m))
											*std::pow(L4,(s-t));
										double factorials =
											Factorial(n)
											*Factorial(m-n)
											*Factorial(t-m)
											*Factorial(s-t);
										double res = orig*powers/factorials;
										
										arr.get()[j][k] += res;
                                    }
								}
							}

							// RecRel #2:
							{
								double fac1 = (
									0.5*(16.0*PI_2*_alpha_s.Beta1() + 4*PI*r1*_alpha_s.Beta2() - (c + b*r1)*_alpha_s.Beta3())
								) * _D2[j][1][s][s][1][k];
								double fac2 = RecRelN3LO_2(_D2[j][0][s-1][s-1][0], k,
                                    getSplitFunc("P0ns"), getSplitFunc("P1nsm"), getSplitFunc("P2nsm"), getSplitFunc("P3nsm"));
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
                                    Factorial(n)
                                    *Factorial(m-n)
                                    *Factorial(t-m)
                                    *Factorial(s-t);
                                double res = orig*powers/factorials;
                                
                                arr.get()[j][k] += res;
							}

							// RecRel #3:
							for (int m=s-1; m>=0; m--)
							{
								double fac1 = -(
									16.0*PI_2*_alpha_s.Beta1() + 4.0*PI*r1*_alpha_s.Beta2() + r1*r1*_alpha_s.Beta3()
								) * _D2[j][1][s][m+1][1][k];
								double fac2 = RecRelN3LO_3(_D2[j][0][s-1][m][0], k,
                                    getSplitFunc("P0ns"), getSplitFunc("P1nsm"), getSplitFunc("P2nsm"), getSplitFunc("P3nsm"));
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
                                    Factorial(n)
                                    *Factorial(m-n)
                                    *Factorial(t-m)
                                    *Factorial(s-t);
                                double res = orig*powers/factorials;
                                
                                arr.get()[j][k] += res;
							}

							// RecRel #4:
							for (int t=s-1; t>=0; t--)
							{
								for (int m=t; m>=0; m--)
								{
									double fac1a = -2.0*b*gamma;
									double fac1b = 32*PI_2*(b+r1)*_alpha_s.Beta1() - 8*PI*c*_alpha_s.Beta2() - 2*c*r1*_alpha_s.Beta3();
									double fac1 = fac1a*_D2[j][1][t+1][m+1][0][k] + fac1b*_D2[j][1][t+1][m+1][1][k];
									double fac2 = RecRelN3LO_4(_D2[j][0][t][m][0], k,
                                        getSplitFunc("P0ns"), getSplitFunc("P1nsm"), getSplitFunc("P2nsm"), getSplitFunc("P3nsm"));
									_D2[j][1][t][m][0][k] = (fac1 + fac2)/gamma;

                                    uint n = 0;
                                    double orig = _D2[j][1][t][m][0][k];
                                    double powers =
                                        std::pow(L1,n)
                                        *std::pow(L2,(m-n))
                                        *std::pow(L3,(t-m))
                                        *std::pow(L4,(s-t));
                                    double factorials =
                                        Factorial(n)
                                        *Factorial(m-n)
                                        *Factorial(t-m)
                                        *Factorial(s-t);
                                    double res = orig*powers/factorials;
                                    
                                    arr.get()[j][k] += res;
								}
							}
						}

						// plus distributions
						for (uint j=32; j<=30+_nf; j++)
						{
						    // RecRel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint m=1; m<=t; m++)
								{
									for (uint n=1; n<=m; n++)
                                    {
										_D2[j][1][t][m][n][k] = RecRelN3LO_1(_D2[j][0][t-1][m-1][n-1], k, getSplitFunc("P0ns"));

                                        double orig = _D2[j][1][t][m][n][k];
										double powers =
											std::pow(L1,n)
											*std::pow(L2,(m-n))
											*std::pow(L3,(t-m))
											*std::pow(L4,(s-t));
										double factorials =
											Factorial(n)
											*Factorial(m-n)
											*Factorial(t-m)
											*Factorial(s-t);
										double res = orig*powers/factorials;
										
										arr.get()[j][k] += res;
                                    }
								}
							}

							// RecRel #2:
							{
								double fac1 = (
									0.5*(16.0*PI_2*_alpha_s.Beta1() + 4*PI*r1*_alpha_s.Beta2() - (c + b*r1)*_alpha_s.Beta3())
								) * _D2[j][1][s][s][1][k];
								double fac2 = RecRelN3LO_2(_D2[j][0][s-1][s-1][0], k,
                                    getSplitFunc("P0ns"), getSplitFunc("P1nsp"), getSplitFunc("P2nsp"), getSplitFunc("P3nsp"));
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
                                    Factorial(n)
                                    *Factorial(m-n)
                                    *Factorial(t-m)
                                    *Factorial(s-t);
                                double res = orig*powers/factorials;
                                
                                arr.get()[j][k] += res;
							}

							// RecRel #3:
							for (int m=s-1; m>=0; m--)
							{
								double fac1 = -(
									16.0*PI_2*_alpha_s.Beta1() + 4.0*PI*r1*_alpha_s.Beta2() + r1*r1*_alpha_s.Beta3()
								) * _D2[j][1][s][m+1][1][k];
								double fac2 = RecRelN3LO_3(_D2[j][0][s-1][m][0], k,
                                    getSplitFunc("P0ns"), getSplitFunc("P1nsp"), getSplitFunc("P2nsp"), getSplitFunc("P3nsp"));
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
                                    Factorial(n)
                                    *Factorial(m-n)
                                    *Factorial(t-m)
                                    *Factorial(s-t);
                                double res = orig*powers/factorials;
                                
                                arr.get()[j][k] += res;
							}

							// RecRel #4:
							for (int t=s-1; t>=0; t--)
							{
								for (int m=t; m>=0; m--)
								{
									double fac1a = -2.0*b*gamma;
									double fac1b = 32*PI_2*(b+r1)*_alpha_s.Beta1() - 8*PI*c*_alpha_s.Beta2() - 2*c*r1*_alpha_s.Beta3();
									double fac1 = fac1a*_D2[j][1][t+1][m+1][0][k] + fac1b*_D2[j][1][t+1][m+1][1][k];
									double fac2 = RecRelN3LO_4(_D2[j][0][t][m][0], k,
                                        getSplitFunc("P0ns"), getSplitFunc("P1nsp"), getSplitFunc("P2nsp"), getSplitFunc("P3nsp"));
									_D2[j][1][t][m][0][k] = (fac1 + fac2)/gamma;

                                    uint n = 0;
                                    double orig = _D2[j][1][t][m][0][k];
                                    double powers =
                                        std::pow(L1,n)
                                        *std::pow(L2,(m-n))
                                        *std::pow(L3,(t-m))
                                        *std::pow(L4,(s-t));
                                    double factorials =
                                        Factorial(n)
                                        *Factorial(m-n)
                                        *Factorial(t-m)
                                        *Factorial(s-t);
                                    double res = orig*powers/factorials;
                                    
                                    arr.get()[j][k] += res;
								}
							}
						}

						// valence distribution
						{
						    // RecRel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint m=1; m<=t; m++)
								{
									for (uint n=1; n<=m; n++)
                                    {
										_D2[25][1][t][m][n][k] = RecRelN3LO_1(_D2[25][0][t-1][m-1][n-1], k, getSplitFunc("P0ns"));

                                        double orig = _D2[25][1][t][m][n][k];
										double powers =
											std::pow(L1,n)
											*std::pow(L2,(m-n))
											*std::pow(L3,(t-m))
											*std::pow(L4,(s-t));
										double factorials =
											Factorial(n)
											*Factorial(m-n)
											*Factorial(t-m)
											*Factorial(s-t);
										double res = orig*powers/factorials;
										
										arr.get()[25][k] += res;
                                    }
								}
							}

							// RecRel #2:
							{
								double fac1 = (
									0.5*(16.0*PI_2*_alpha_s.Beta1() + 4*PI*r1*_alpha_s.Beta2() - (c + b*r1)*_alpha_s.Beta3())
								) * _D2[25][1][s][s][1][k];
								double fac2 = RecRelN3LO_2(_D2[25][0][s-1][s-1][0], k,
                                    getSplitFunc("P0ns"), getSplitFunc("P1nsm"), getSplitFunc("P2nsv"), getSplitFunc("P3nsv"));
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
                                    Factorial(n)
                                    *Factorial(m-n)
                                    *Factorial(t-m)
                                    *Factorial(s-t);
                                double res = orig*powers/factorials;
                                
                                arr.get()[25][k] += res;
							}

							// RecRel #3:
							for (int m=s-1; m>=0; m--)
							{
								double fac1 = -(
									16.0*PI_2*_alpha_s.Beta1() + 4.0*PI*r1*_alpha_s.Beta2() + r1*r1*_alpha_s.Beta3()
								) * _D2[25][1][s][m+1][1][k];
								double fac2 = RecRelN3LO_3(_D2[25][0][s-1][m][0], k,
                                    getSplitFunc("P0ns"), getSplitFunc("P1nsm"), getSplitFunc("P2nsv"), getSplitFunc("P3nsv"));
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
                                    Factorial(n)
                                    *Factorial(m-n)
                                    *Factorial(t-m)
                                    *Factorial(s-t);
                                double res = orig*powers/factorials;
                                
                                arr.get()[25][k] += res;
							}

							// RecRel #4:
							for (int t=s-1; t>=0; t--)
							{
								for (int m=t; m>=0; m--)
								{
									double fac1a = -2.0*b*gamma;
									double fac1b = 32*PI_2*(b+r1)*_alpha_s.Beta1() - 8*PI*c*_alpha_s.Beta2() - 2*c*r1*_alpha_s.Beta3();
									double fac1 = fac1a*_D2[25][1][t+1][m+1][0][k] + fac1b*_D2[25][1][t+1][m+1][1][k];
									double fac2 = RecRelN3LO_4(_D2[25][0][t][m][0], k,
                                        getSplitFunc("P0ns"), getSplitFunc("P1nsm"), getSplitFunc("P2nsv"), getSplitFunc("P3nsv"));
									_D2[25][1][t][m][0][k] = (fac1 + fac2)/gamma;

                                    uint n = 0;
                                    double orig = _D2[25][1][t][m][0][k];
                                    double powers =
                                        std::pow(L1,n)
                                        *std::pow(L2,(m-n))
                                        *std::pow(L3,(t-m))
                                        *std::pow(L4,(s-t));
                                    double factorials =
                                        Factorial(n)
                                        *Factorial(m-n)
                                        *Factorial(t-m)
                                        *Factorial(s-t);
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



    void DGLAPSolver::FixDistributions2(bool resum_tab, bool resum_threshold, std::vector<ArrayGrid>& temp_arr, std::vector<ArrayGrid>& temp_arr_singlet)
    {
        for (uint t=1; t<=_trunc_idx; ++t)
        {
            for (uint j=0; j<=1; ++j)
            {
                for (uint n=0; n<=1; ++n)
                    _S2[t][j][n].zero();
            }
        }

        double Nf = static_cast<double>(_nf);
        if (resum_tab)
        {
            for (uint k=0; k<_grid.Size()-1;k++)
            {
                if (_order>=2)
                {
                    _F2[13][k]=_F2[25][k];
                    for (uint j=26; j<=24+_nf; j++)
                        _F2[13][k] += _F2[j][k];
                    _F2[13][k] /= Nf;
                    for (uint j=14; j<=12+_nf; j++)
                        _F2[j][k] = _F2[13][k] - _F2[j+12][k];
                }

                _F2[19][k] = _F2[31][k];
                for (uint j=32; j<=30+_nf; j++)
                    _F2[19][k] += _F2[j][k];
                _F2[19][k] /= Nf;

                for (uint j=20; j<=18+_nf; j++)
                    _F2[j][k] = _F2[19][k] - _F2[j+12][k];

                for (uint j=1; j<=_nf; j++)
                {
                    _F2[j][k]   =0.5*(_F2[j+18][k] + _F2[j+12][k]);
                    _F2[j+6][k] =0.5*(_F2[j+18][k] - _F2[j+12][k]);
                }

                if (_order<2)
                {
                    _F2[25][k]=0.0;
                    for (uint j=13; j<=12+_nf; j++)
                        _F2[25][k] += _F2[j][k];

                    for (uint j=26; j<=24+_nf; j++)
                        _F2[j][k] = _F2[13][k] - _F2[j-12][k];
                }
            }
        }
        else if (resum_threshold)
        {
            switch (_order)
            {
                case 0:
                case 1:
                {
                    // if we resummed to a threshold energy,
                    // we must fix the temporary array
                    for (uint k=0; k<_grid.Size()-1; k++)
                    {
                        temp_arr[19][k] = temp_arr_singlet[31][k];
                        for (uint j=32; j<=30+_nf; j++)
                            temp_arr[19][k] += temp_arr[j][k];
                        temp_arr[19][k] /= Nf;

                        for (uint j=20; j<=18+_nf; j++)
                            temp_arr[j][k] = temp_arr[19][k] - temp_arr[j+12][k];

                        for (uint j=1; j<=_nf; j++)
                        {
                            temp_arr[j][k]   = 0.5*(temp_arr[j+18][k] + temp_arr[j+12][k]);
                            temp_arr[j+6][k] = 0.5*(temp_arr[j+18][k] - temp_arr[j+12][k]);
                        }
                    }
                }; break;
                case 2:
                case 3:
                {
                    for (uint k=0; k<_grid.Size()-1;k++)
                    {
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

                        for (uint j=1; j<=_nf; j++)
                        {
                            temp_arr[j][k]  =0.5*(temp_arr[j+18][k]+temp_arr[j+12][k]);
                            temp_arr[j+6][k]=0.5*(temp_arr[j+18][k]-temp_arr[j+12][k]);
                        }
                    }
                }
            }
        }
    }

    void DGLAPSolver::HeavyFlavorTreatment2()
    {
        std::cerr << "[DGLAP2] Treating heavy flavors: "
				  << _nf+1
				  << "th quark mass threshold (mass "
				  << _alpha_s.Masses(_nf+1) << ")\n";

		OpMatElem::update(-_log_mur2_muf2, _nf);

		// Copy of pre-threshold distributions
		// the nf+1 dists are defined in terms of the nf dists,
		// so we need this copy since we would otherwise be overwriting
		// a given dist as we do the calculations
		// we just store them in the s=1 array, and modify the s=0 array
		// (which are the initial conditions for the next set of iterations
		// at the next nf
		std::cerr << "[DGLAP] Creating copy of pre-threshold distributions... ";
        
        std::vector<ArrayGrid> arr(13, ArrayGrid{_grid});
        std::vector<ArrayGrid> arr_singlet(2, ArrayGrid{_grid});

		for (uint j=0; j<=1; ++j)
			arr_singlet[j] = _S2[0][j][0];
		
		if (_order == 2) 
		{
			for (uint i=1; i<=_nf; i++)
			{
				for (uint j=i; j<=i+6; j+=6)
					arr[j] = _C2[j][0][0][0];
			}
		} else if (_order == 3) {
            for (uint i=1; i<=_nf; i++)
			{
				for (uint j=i; j<=i+6; j+=6)
					arr[j] = _D2[j][0][0][0][0];
			}
        }
		std::cerr << "Done.\n";

		double as = _alpha_s.Post(_nf+1);
		std::cerr << "[DGLAP2] Value of alpha_s post threshold: " << as << '\n';

		switch(_order)
		{
			case 2:
			{
				for (uint k=0; k<_grid.Size()-1;k++)
				{
					// q
					for (uint j=1; j<=_nf; j++)
						HFT2_NNLO1(arr[j], j, k);
					// qbar
					for (uint j=1+6; j<=_nf+6; j++)
						HFT2_NNLO1(arr[j], j, k);

					HFT2_NNLO2(arr_singlet[1], arr_singlet[0], k); // gluon
					HFT2_NNLO3(arr_singlet[1], arr_singlet[0], k); // heavy flavor
				}
			} break;
            case 3:
            {
                for (uint k=0; k<_grid.Size()-1;k++)
				{
                    const double fac_n3lo = as*as*as/(64.0*PI_3);
					const double convSPa = getSplitFunc("A3psqq").convolution(arr_singlet[1], k);
					const double convSPb = getSplitFunc("A3sqg").convolution(arr_singlet[0], k);
					const double SP = fac_n3lo*(convSPa + convSPb)/static_cast<double>(_nf);

					// q
					for (uint j=1; j<=_nf; j++)
						HFT2_N3LO1(arr[j], arr[j+6], j, k, SP);
					// qbar
					for (uint j=1+6; j<=_nf+6; j++)
						HFT2_N3LO2(arr[j-6], arr[j], j, k, SP);

					HFT2_N3LO3(arr_singlet[0], arr_singlet[1], k); // gluon
					HFT2_N3LO4(arr_singlet[0], arr_singlet[1], k); // heavy flavor
				}
            } break;
		}
    }
    void DGLAPSolver::HFT2_NNLO1(ArrayGrid& c, uint j, uint k)
    {
        double const as = _alpha_s.Post(_nf+1);
        double const conv = getSplitFunc("A2ns").convolution(c, k);
        
		_C2[j][0][0][0][k] += std::pow(as/(4.0*PI), 2) * conv;
    }
    void DGLAPSolver::HFT2_NNLO2(ArrayGrid& s1, ArrayGrid& s2, uint k)
    {
        double const as = _alpha_s.Post(_nf+1);
        double const conv1 = getSplitFunc("A2gq").convolution(s1, k);
        double const conv2 = getSplitFunc("A2gg").convolution(s2, k);

		_S2[0][0][0][k] += std::pow(as/(4.0*PI), 2) * (conv1 + conv2);
    }
    void DGLAPSolver::HFT2_NNLO3(ArrayGrid& s1, ArrayGrid& s2, uint k)
    {
        double const as = _alpha_s.Post(_nf+1);
        double const conv1 = getSplitFunc("A2hq").convolution(s1, k);
        double const conv2 = getSplitFunc("A2hg").convolution(s2, k);
		double const fac = 0.5*std::pow(as/(4.0*PI), 2) * (conv1 + conv2);

		_C2[_nf+1][0][0][0][k] = fac;
        _C2[_nf+7][0][0][0][k] = fac;
    }

    // q
	void DGLAPSolver::HFT2_N3LO1(ArrayGrid& q, ArrayGrid& qb, uint j, uint k, double SP)
	{
	    const double as = _alpha_s.Post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

		/*
		const double conv1a = getSplitFunc("A2ns").convolution(q, k);
		const double conv1b = getSplitFunc("A3nsp").convolution(q, k);
		const double conv1c = conv1a;
		const double conv1d = getSplitFunc("A3nsm").convolution(q, k);

		const double conv2a = getSplitFunc("A2ns").convolution(qb, k);
		const double conv2b = getSplitFunc("A3nsp").convolution(qb, k);
		const double conv2c = conv2a;
		const double conv2d = getSplitFunc("A3nsm").convolution(qb, k);
		
		_D2[j][0][0][0][0][k] += 0.5*(
			  ((fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv1c + fac_n3lo*conv1d))
			+ ((fac_nnlo*conv2a + fac_n3lo*conv2b) - (fac_nnlo*conv2c + fac_n3lo*conv2d))
			+ SP);
		*/
		
		double const conv1a = getSplitFunc("A2ns").convolution(q, k);
		double const conv1b = getSplitFunc("A3nsp").convolution(q, k, true);

		double const conv2a = getSplitFunc("A2ns").convolution(qb, k);
		double const conv2b = getSplitFunc("A3nsp").convolution(qb, k, true);

		double const conv3a = conv1a;
		double const conv3b = getSplitFunc("A3nsm").convolution(q, k, true);

		double const conv4a = conv2a;
		double const conv4b = getSplitFunc("A3nsm").convolution(qb, k, true);

		_D2[j][0][0][0][0][k] += 0.5*(
			((fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b)) +
			((fac_nnlo*conv3a + fac_n3lo*conv3b) - (fac_nnlo*conv4a + fac_n3lo*conv4b)) +
			SP);
	}

    // qbar
	void DGLAPSolver::HFT2_N3LO2(ArrayGrid& q, ArrayGrid& qb, uint j, uint k, double SP)
	{
		const double as = _alpha_s.Post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

		/*
		const double conv1a = getSplitFunc("A2ns").convolution(q, k);
		const double conv1b = getSplitFunc("A3nsp").convolution(q, k);
		const double conv1c = conv1a;
		const double conv1d = getSplitFunc("A3nsm").convolution(q, k);

		const double conv2a = getSplitFunc("A2ns").convolution(qb, k);
		const double conv2b = getSplitFunc("A3nsp").convolution(qb, k);
		const double conv2c = conv2a;
		const double conv2d = getSplitFunc("A3nsm").convolution(qb, k);

		_D2[j][0][0][0][0][k] += 0.5*(
			((fac_nnlo*conv1a + fac_n3lo*conv1b) - (fac_nnlo*conv1c + fac_n3lo*conv1d))
			+ ((fac_nnlo*conv2a + fac_n3lo*conv2b) + (fac_nnlo*conv2c + fac_n3lo*conv2d))
			+ SP);
		*/

		double const conv1a = getSplitFunc("A2ns").convolution(q, k);
		double const conv1b = getSplitFunc("A3nsp").convolution(q, k, true);

		double const conv2a = getSplitFunc("A2ns").convolution(qb, k);
		double const conv2b = getSplitFunc("A3nsp").convolution(qb, k, true);

		double const conv3a = conv1a;
		double const conv3b = getSplitFunc("A3nsm").convolution(q, k, true);

		double const conv4a = conv2a;
		double const conv4b = getSplitFunc("A3nsm").convolution(qb, k, true);

		_D2[j][0][0][0][0][k] += 0.5*(
			((fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b)) -
			((fac_nnlo*conv3a + fac_n3lo*conv3b) - (fac_nnlo*conv4a + fac_n3lo*conv4b)) +
			SP);
	}

	// gluon (index 0 in S array)
	void DGLAPSolver::HFT2_N3LO3(ArrayGrid& g, ArrayGrid& qp, uint k)
	{
		const double as = _alpha_s.Post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);
	    
		const double conv1a = getSplitFunc("A2gq").convolution(qp, k);
		const double conv1b = getSplitFunc("A3gq").convolution(qp, k, true);
		const double conv2a = getSplitFunc("A2gg").convolution(g, k);
		const double conv2b = getSplitFunc("A3gg").convolution(g, k, true);
		
		_S2[0][0][0][k] += (fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b);
	}

	// heavy flavor
	void DGLAPSolver::HFT2_N3LO4(ArrayGrid& g, ArrayGrid& qp, uint k)
	{
		const double as = _alpha_s.Post(_nf+1);
		const double fac_nnlo = as*as/(16.0*PI_2);
		const double fac_n3lo = as*as*as/(64.0*PI_3);

	    const double conv1a = getSplitFunc("A2hq").convolution(qp, k);
		const double conv1b = getSplitFunc("A3hq").convolution(qp, k, true);
		const double conv2a = getSplitFunc("A2hg").convolution(g, k);
		const double conv2b = getSplitFunc("A3hg").convolution(g, k, true);
		
        const double res = 0.5*((fac_nnlo*conv1a + fac_n3lo*conv1b) + (fac_nnlo*conv2a + fac_n3lo*conv2b));
		_D2[(_nf+1)][0][0][0][0][k] = res;
        _D2[(_nf+1)+6][0][0][0][0][k] = res;
	}






    void DGLAPSolver::EvolveNonSinglet2Threaded(
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			double L1, double L2, double L3, double L4)
	{
		switch (_order)
		{
			case 0: // LO
			{
				std::cout << "[THREAD2] Performing LO non-singlet evolution threaded.\n";

                for (uint j=13; j<=12+_nf; ++j)
                    arr.get()[j] = _A2[j][0];
                for (uint j=32; j<=30+_nf; ++j)
                    arr.get()[j] = _A2[j][0];
                    
				std::vector<std::thread> threads{};

				for (uint j=13; j<=12+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt2_EvolveDistribution_NS_LO, this, arr, j, L1);	
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt2_EvolveDistribution_NS_LO, this, arr, j, L1);

			    for (std::thread & t : threads)
					t.join();

				std::cout << "[THREAD2] Finished performing threaded LO non-singlet evolution.\n";
			} break;
			case 1: // NLO
			{
				std::cout << "[THREAD2] Performing NLO non-singlet evolution threaded.\n";

                for (uint j=13; j<=12+_nf; ++j)
                    arr.get()[j] = _B2[j][0][0];
                for (uint j=32; j<=30+_nf; ++j)
                    arr.get()[j] = _B2[j][0][0];

				std::vector<std::thread> threads{};
                std::array<double, 2> L{L1, L2};

				for (uint j=13; j<=12+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt2_EvolveDistribution_NS_NLO, this, arr, j, "P1nsm", L);
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt2_EvolveDistribution_NS_NLO, this, arr, j, "P1nsp", L);

				for (std::thread & t : threads)
					t.join();

				std::cout << "[THREAD2] Finished performing threaded NLO non-singlet evolution.\n";
			} break;
			case 2: // NNLO
			{
				std::cout << "[THREAD2] Performing NNLO non-singlet evolution threaded.\n";
				std::vector<std::thread> threads{};

                std::array<double, 3> L{L1, L2, L3};
				std::array<std::string, 2> nsm{"P1nsm", "P2nsm"};
				std::array<std::string, 2> nsp{"P1nsp", "P2nsp"};
				std::array<std::string, 2> nsv{"P1nsm", "P2nsv"};

				for (uint j=26; j<=24+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt2_EvolveDistribution_NS_NNLO, this, arr, j, nsm, L);
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt2_EvolveDistribution_NS_NNLO, this, arr, j, nsp, L);
				threads.emplace_back(&DGLAPSolver::_mt2_EvolveDistribution_NS_NNLO, this, arr, 25, nsv, L);
				
				for (std::thread & t : threads)
					t.join();

				std::cout << "[THREAD2] Finished performing threaded NNLO non-singlet evolution.\n";
			} break;
			case 3: // N3LO nonsinglet
			{
				std::cout << "[THREAD2] Performing N3LO non-singlet evolution threaded.\n";

				for (uint j=26; j<=24+_nf; ++j)
                    arr.get()[j] = _D2[j][0][0][0][0];
                for (uint j=32; j<=30+_nf; ++j)
                    arr.get()[j] = _D2[j][0][0][0][0];
				arr.get()[25] = _D2[25][0][0][0][0];
				
				std::vector<std::thread> threads{};

                std::array<double, 4> L{L1, L2, L3, L4};
				std::array<std::string, 3> nsm{"P1nsm", "P2nsm", "P3nsm"};
				std::array<std::string, 3> nsp{"P1nsp", "P2nsp", "P3nsp"};
				std::array<std::string, 3> nsv{"P1nsm", "P2nsv", "P3nsv"};

				std::cout << "=============== BEGIN THREAD OUTPUT ====================\n";
				for (uint j=26; j<=24+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt2_EvolveDistribution_NS_N3LO, this, arr, j, nsm, L);
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt2_EvolveDistribution_NS_N3LO, this, arr, j, nsp, L);
				threads.emplace_back(&DGLAPSolver::_mt2_EvolveDistribution_NS_N3LO, this, arr, 25, nsv, L);
				
				for (std::thread & t : threads)
					t.join();

				std::cout << "=============== END THREAD OUTPUT ====================\n";
				std::cout << "[THREAD2] Finished performing threaded N3LO non-singlet evolution.\n";
			} break;
		}
	}


    void DGLAPSolver::_mt2_EvolveDistribution_NS_LO (
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
        uint j, double L1)
    {
        for (uint n=0; n<_iterations-1; n++)
		{
			for (uint k=0; k<_grid.Size()-1; k++)
            {
				_A[j][1][k] = RecRelLO(_A[j][0], k, _P0ns);
                arr.get()[j][k] += _A2[j][1][k]*std::pow(L1, n)/Factorial(n);
            }
			_A2[j][0] = _A2[j][1];
		}
		std::cout << std::format("[THREAD2] LO Distribution j={} finished.\n", j);
    }
    void DGLAPSolver::_mt2_EvolveDistribution_NS_NLO(
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::string const& P1, std::array<double, 2> const&  L)
    {
        double const L1 = L[0];
        double const L2 = L[1];
        for (uint s=1; s<_iterations; s++)
        {
            for (uint k=0; k<_grid.Size()-1;k++)
            {
                for (uint n=1; n<=s; n++)
                {
                    _B2[j][1][n][k] = RecRelNLO_1(_B2[j][0][n-1], k, getSplitFunc("P0ns"));
                    arr.get()[j][k] += _B2[j][1][n][k]*std::pow(L1,n)*std::pow(L2,s-n)/Factorial(n)/Factorial(s-n);
                }

                uint n = 0;
                double res = RecRelNLO_2(_B2[j][0][0], k, 
                    getSplitFunc("P0ns"), 
                    getSplitFunc(P1));
                _B2[j][1][0][k] = -_B2[j][1][1][k] + res;
                arr.get()[j][k] += _B2[j][1][0][k]
                    *std::pow(L1,n)*std::pow(L2,s-n)
                    /Factorial(n)/Factorial(s-n);
            }
            for (uint n=0; n<=s; ++n)
                _B2[j][0][n] = _B2[j][1][n];
        }
    }
    void DGLAPSolver::_mt2_EvolveDistribution_NS_NNLO(
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::array<std::string, 2> const& P, std::array<double, 3> const& L)
    {
        double const L1 = L[0];
        double const L2 = L[1];
        double const L3 = L[2];

        for (uint s=1; s<_iterations; s++)
        {
            for (uint k=0; k<_grid.Size()-1; k++)
            {
                // RecRel #1:
                for (uint t=1; t<=s; t++)
                {
                    for (uint n=1; n<=t; n++)
                    {
                        double recrel = RecRelNNLO_1(_C2[j][0][t-1][n-1], k, getSplitFunc("P0ns"));
                        _C2[j][1][t][n][k] = recrel;

                        double orig = _C2[j][1][t][n][k];
                        double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                        double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
                        double res = orig*powers/factorials;

                        arr.get()[j][k] += res;
                    }
                }

                // RecRel #2:
                {
                    double fac1 = -0.5*_C2[j][1][s][1][k];
                    double fac2 = RecRelNNLO_2(_C2[j][0][s-1][0], k, 
                        getSplitFunc("P0ns"), getSplitFunc(P[0]), getSplitFunc(P[1]));
                    _C2[j][1][s][0][k] = fac1 + fac2;

                    uint n = 0;
                    uint t = s;
                    double orig = _C2[j][1][s][0][k];
                    double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                    double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
                    double res = orig*powers/factorials;

                    arr.get()[j][k] += res;
                }

                // these must be regular ints;
                // unsigned ints, when they are 0 and get --,
                // underflow back to positive 4b,
                // remaining positive and the loop continues (very bad!)

                // RecRel #3:
                for (int t=s-1; t>=0; t--)
                {
                    double fac1 = -2.0*_alpha_s.Beta1()*(_C2[j][1][t+1][0][k] + _C2[j][1][t+1][1][k]);
                    double fac2 = RecRelNNLO_3(_C2[j][0][t][0], k, 
                        getSplitFunc("P0ns"), getSplitFunc(P[0]));
                    _C2[j][1][t][0][k] = fac1 + fac2;

                    uint n = 0;
                    double orig = _C2[j][1][t][0][k];
                    double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                    double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
                    double res = orig*powers/factorials;

                    arr.get()[j][k] += res;
                }
            }

            for (uint t=0; t<=s; ++t)
            {
                for (uint n=0; n<=t; ++n)
                    _C2[j][0][t][n] = _C2[j][1][t][n];
            }
        }
    }
    void DGLAPSolver::_mt2_EvolveDistribution_NS_N3LO(
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::array<std::string, 3> const& P, std::array<double, 4> const& L)
    {
		std::println("  [j={}] Beginning evolution...", j);
		
        double const L1 = L[0];
        double const L2 = L[1];
        double const L3 = L[2];
        double const L4 = L[3];

        // some shorthand
        double r1 = _r1[_nf];
        double b = _b[_nf];
        double c = _c[_nf];
        double gamma = (r1*r1 + r1*b + c)*_alpha_s.Beta3();

        for (uint s=1; s<_iterations; s++)
        {
			std::println("  [j={}] Iteration {}/{}", j, s, _iterations-1);
            for (uint k=0; k<_grid.Size()-1; k++)
            {
                // RecRel #1:
                for (uint t=1; t<=s; t++)
                {
                    for (uint m=1; m<=t; m++)
                    {
                        for (uint n=1; n<=m; n++)
                        {
                            _D2[j][1][t][m][n][k] = RecRelN3LO_1(_D2[j][0][t-1][m-1][n-1], k, getSplitFunc("P0ns"));

                            double orig = _D2[j][1][t][m][n][k];
                            double powers =
                                std::pow(L1,n)
                                *std::pow(L2,(m-n))
                                *std::pow(L3,(t-m))
                                *std::pow(L4,(s-t));
                            double factorials =
                                Factorial(n)
                                *Factorial(m-n)
                                *Factorial(t-m)
                                *Factorial(s-t);
                            double res = orig*powers/factorials;
                            
                            arr.get()[j][k] += res;
                        }
                    }
                }

                // RecRel #2:
                {
                    double fac1 = (
                        0.5*(16.0*PI_2*_alpha_s.Beta1() + 4*PI*r1*_alpha_s.Beta2() - (c + b*r1)*_alpha_s.Beta3())
                    ) * _D2[j][1][s][s][1][k];
                    double fac2 = RecRelN3LO_2(_D2[j][0][s-1][s-1][0], k,
                        getSplitFunc("P0ns"), getSplitFunc(P[0]), getSplitFunc(P[1]), getSplitFunc(P[2]));
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
                        Factorial(n)
                        *Factorial(m-n)
                        *Factorial(t-m)
                        *Factorial(s-t);
                    double res = orig*powers/factorials;
                    
                    arr.get()[j][k] += res;
                }

                // RecRel #3:
                for (int m=s-1; m>=0; m--)
                {
                    double fac1 = -(
                        16.0*PI_2*_alpha_s.Beta1() + 4.0*PI*r1*_alpha_s.Beta2() + r1*r1*_alpha_s.Beta3()
                    ) * _D2[j][1][s][m+1][1][k];
                    double fac2 = RecRelN3LO_3(_D2[j][0][s-1][m][0], k,
                        getSplitFunc("P0ns"), getSplitFunc(P[0]), getSplitFunc(P[1]), getSplitFunc(P[2]));
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
                        Factorial(n)
                        *Factorial(m-n)
                        *Factorial(t-m)
                        *Factorial(s-t);
                    double res = orig*powers/factorials;
                    
                    arr.get()[j][k] += res;
                }

                // RecRel #4:
                for (int t=s-1; t>=0; t--)
                {
                    for (int m=t; m>=0; m--)
                    {
                        double fac1a = -2.0*b*gamma;
                        double fac1b = 32*PI_2*(b+r1)*_alpha_s.Beta1() - 8*PI*c*_alpha_s.Beta2() - 2*c*r1*_alpha_s.Beta3();
                        double fac1 = fac1a*_D2[j][1][t+1][m+1][0][k] + fac1b*_D2[j][1][t+1][m+1][1][k];
                        double fac2 = RecRelN3LO_4(_D2[j][0][t][m][0], k,
                            getSplitFunc("P0ns"), getSplitFunc(P[0]), getSplitFunc(P[1]), getSplitFunc(P[2]));
                        _D2[j][1][t][m][0][k] = (fac1 + fac2)/gamma;

                        uint n = 0;
                        double orig = _D2[j][1][t][m][0][k];
                        double powers =
                            std::pow(L1,n)
                            *std::pow(L2,(m-n))
                            *std::pow(L3,(t-m))
                            *std::pow(L4,(s-t));
                        double factorials =
                            Factorial(n)
                            *Factorial(m-n)
                            *Factorial(t-m)
                            *Factorial(s-t);
                        double res = orig*powers/factorials;
                        
                        arr.get()[j][k] += res;
                    }
                }
            }

            for (uint t=0; t<=s; ++t)
            {
                for (uint m=0; m<=s; ++m)
                {
                    for (uint n=0; n<=s; ++n)
                        _D2[j][0][t][m][n] = _D2[j][1][t][m][n];
                }
            }
        }
    }
}
