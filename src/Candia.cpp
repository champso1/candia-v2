#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/FuncArrayGrid.hpp"
#include "Candia-v2/Grid.hpp"
#include "Candia-v2/SplittingFn.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"

#include <cstdlib>
#include <functional>
#include <print>



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
		uint order, Grid & grid, double Qf,
		uint iterations, uint trunc_idx,
		std::unique_ptr<Distribution> initial_dist,
		double kr, bool multi_thread) 
		: _order{order},  _grid{grid}, _Qf{Qf},
		  _alpha_s{order, initial_dist->Q0(), initial_dist->alpha0(), initial_dist->masses(), kr},
		  _mur2_muf2{kr}, _log_mur2_muf2{std::log(kr)},
		  _dist{std::move(initial_dist)},
		  _iterations{iterations}, _trunc_idx{trunc_idx},
		  _multi_thread{multi_thread}
	{
		std::println("[DGLAP] Evolving with log(mu_R / mu_F) = log({:.1}) = {:.4}.", _mur2_muf2, _log_mur2_muf2);
		
		std::print("[DGLAP] Reserving space in coefficient arrays...");
		switch(_order)
		{
			case 0:
			{
				_trunc_idx = 0; // LO has exact singlet solution, do not add additional terms

				_A2 = std::vector<std::vector<ArrayGrid>>{
					DISTS, std::vector<ArrayGrid>{
						2, ArrayGrid{grid}
					}
				};
			} break;
			case 1:
			{
				_B2 = MultiDimArrayGrid_t<3>{
					DISTS, MultiDimArrayGrid_t<2>{
						2, MultiDimArrayGrid_t<1>{
							_iterations, ArrayGrid{grid}
						}
					}
				};

			} break;
			case 2:
			{
				_C2 = MultiDimArrayGrid_t<4>{
					DISTS, MultiDimArrayGrid_t<3>{
						2, MultiDimArrayGrid_t<2>{
							_iterations, MultiDimArrayGrid_t<1>{
								_iterations, ArrayGrid{grid}
							}
						}
					}
				};
			} break;
			case 3:
			{
				_D2 = MultiDimArrayGrid_t<5>{
					DISTS, MultiDimArrayGrid_t<4>{
						2, MultiDimArrayGrid_t<3>{
							_iterations, MultiDimArrayGrid_t<2>{
								_iterations, MultiDimArrayGrid_t<1>{
									_iterations, ArrayGrid{grid}
								}
							}
						}
					}
				};
			
				// also initialize the N3LO coefficients we will need
				_r1[3] = -1.11203;
				_b[3] = -2.0 * 0.221422;
				_c[3] = std::pow(0.221422, 2) + std::pow(1.13108, 2);

				_r1[4] = -1.20902;
				_b[4] = -2.0 * 0.286759;
				_c[4] = std::pow(0.286759, 2) + std::pow(1.27279, 2);
				
				_r1[5] = -1.32059;
				_b[5] = -2.0 * 0.42477;
				_c[5] = std::pow(0.42477, 2) + std::pow(1.48548, 2);
				
				_r1[6] = -1.4278;
				_b[6] = -2.0 * 0.796497;
				_c[6] = std::pow(0.796497, 2) + std::pow(1.81681, 2);
			} break;
			default:
			{
				std::println("\n[DGLAP: ERROR] DGLAPSolver(): Found {} for the order, expected a value in range [0, 3].", order);
				exit(EXIT_FAILURE);
			}
		}
		_S2 = decltype(_S2){
			trunc_idx+1, std::vector<std::vector<ArrayGrid>>{
				2, std::vector<ArrayGrid>{
					2, ArrayGrid{_grid}
				}
			}
		};
		std::println(" Done.");


		_F2 = std::vector<ArrayGrid>{
			DISTS, ArrayGrid{grid}
		};

		_alpha_s.calculateThresholdValues(Qf);
		setInitialConditions();
	    _nfi = _dist->nfi();
		_nff = _alpha_s.nff(_nfi, _Qf);
		std::println("[DGLAP] Evolving to {} flavors.", _nff);

		_qct = 0;
	}

	DGLAPSolver::~DGLAPSolver()
	{
		std::println("[DGLAP] Exiting...");
	}


	FunctionGrid& DGLAPSolver::getSplitFunc(std::string const& name)
	{
		auto splitfunc = _expression_grids.find(name);
		if (splitfunc == _expression_grids.end())
		{
			std::println("[DGLAP] Splitting function '{}' does not exist.", name);
			exit(EXIT_FAILURE);
		}
		return _expression_grids.at(name);
	}

	FunctionGrid& DGLAPSolver::getOME(std::string const& name)
	{
		return _expression_grids.at(name);
	}

	void DGLAPSolver::setInitialConditions()
	{
		std::print("[DGLAP] Setting initial conditions... ");
	   
		for (uint k=0; k<_grid.size()-1; k++)
		{
			double x = _grid[k];
			_S2[0][0][0][k] = _dist->xg(x);
			_S2[0][1][0][k] =
				_dist->xuv(x)
				+ 2.0*_dist->xub(x)
				+ _dist->xdv(x)
				+ 2.0*_dist->xdb(x)
				+ 2.0*_dist->xs(x);
		}
	    
		switch (_order)
		{
			case 0:
			{
				for (uint k=0; k<_grid.size()-1; k++)
				{
					double x = _grid[k];
					_A2[7][0][k] = _dist->xub(x);
					_A2[1][0][k] = _dist->xuv(x) + _A2[7][0][k];
					_A2[8][0][k] = _dist->xdb(x);
					_A2[2][0][k] = _dist->xdv(x) + _A2[8][0][k];
					_A2[3][0][k] = _dist->xs(x);
					_A2[9][0][k] = _dist->xs(x);
				}
			} break;
			case 1:
			{
				for (uint k=0; k<_grid.size()-1; k++)
				{
					double x = _grid[k];
					_B2[7][0][0][k] = _dist->xub(x);
					_B2[1][0][0][k] = _dist->xuv(x) + _B2[7][0][0][k];
					_B2[8][0][0][k] = _dist->xdb(x);
					_B2[2][0][0][k] = _dist->xdv(x) + _B2[8][0][0][k];
					_B2[3][0][0][k] = _dist->xs(x);
					_B2[9][0][0][k] = _dist->xs(x);
				}
			} break;
			case 2:
			{
				for (uint k=0; k<_grid.size()-1; k++)
				{
					double x = _grid[k];
					_C2[7][0][0][0][k] = _dist->xub(x);
					_C2[1][0][0][0][k] = _dist->xuv(x) + _C2[7][0][0][0][k];
					_C2[8][0][0][0][k] = _dist->xdb(x);
					_C2[2][0][0][0][k] = _dist->xdv(x) + _C2[8][0][0][0][k];
					_C2[3][0][0][0][k] = _dist->xs(x);
					_C2[9][0][0][0][k] = _dist->xs(x);
				}
			} break;
			case 3:
			{
				for (uint k=0; k<_grid.size()-1; k++)
				{
					double x = _grid[k];
					_D2[7][0][0][0][0][k] = _dist->xub(x);
					_D2[1][0][0][0][0][k] = _dist->xuv(x) + _D2[7][0][0][0][0][k];
					_D2[8][0][0][0][0][k] = _dist->xdb(x);
					_D2[2][0][0][0][0][k] = _dist->xdv(x) + _D2[8][0][0][0][0][k];
					_D2[3][0][0][0][0][k] = _dist->xs(x);
					_D2[9][0][0][0][0][k] = _dist->xs(x);
				}

			} break;
		}

		std::println("Done.");
	}

	void DGLAPSolver::setEvolutionVariables(uint iterations, uint trunc_idx)
	{
		if (iterations > 15)
			std::println("[DGLAP: WARNING] SetEvolutionVariables(): {} iterations is large.", iterations);
		if (trunc_idx > 30)
			std::println("[DGLAP: WARNING] SetEvolutionVariables(): {} truncation terms is large.", trunc_idx);

		_iterations = iterations;
		_trunc_idx = trunc_idx;
	}

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

    void DGLAPSolver::setupCoefficients()
    {
        switch (_order)
		{
			case 0: // LO
			{
				for (uint k=0; k<_grid.size(); k++)
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
				for (uint k=0; k<_grid.size(); k++)
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
				for (uint k=0; k<_grid.size(); k++)
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
				for (uint k=0; k<_grid.size(); k++)
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

	void DGLAPSolver::fixDistributions(bool resum_tab, bool resum_threshold, std::vector<ArrayGrid>& temp_arr, std::vector<ArrayGrid>& temp_arr_singlet)
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
            for (uint k=0; k<_grid.size()-1;k++)
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
                    for (uint k=0; k<_grid.size()-1; k++)
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
                    for (uint k=0; k<_grid.size()-1;k++)
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

	auto DGLAPSolver::evolve() -> decltype(_F2)
	{
		using out_type = decltype(_F2);
		loadAllExpressions();

		std::array<double,1> Qtab{_Qf};
		out_type final_dists;

		// temp array for the threshold summation
		// originally, we wrote directly to the n=0 piece
		// of the coefficient array,
		// but now we only have the two pieces that continually
		// evolve forward
		// so we stick everything into this temp array
		// during the evolution, then move it into the n=0
		// part after the full evolution
		std::vector<ArrayGrid> temp_arr(DISTS, ArrayGrid{_grid});
		std::vector<ArrayGrid> temp_arr_singlet(DISTS, ArrayGrid{_grid});
			
		// since the only difference during the evolution/resummation to
		// the tabulated energy or the threshold is what array we append to, 
		// I create two reference arrays (one for singlet one for non-singlet)
		// that are updated depending on what we are evolving to
		// that way we don't have to evaluate if's inside the grid loop
		// (and just is a bit more readable anyway, assuming I choose a good variable name)
		std::reference_wrapper<std::vector<ArrayGrid>> arr{std::ref(temp_arr)};
		std::reference_wrapper<std::vector<ArrayGrid>> arr_singlet{std::ref(temp_arr_singlet)};

		for (_nf=_nfi; ; _nf++)
		{
			std::println("[DGLAP] Setting nf={}", _nf);

			std::println("[DGLAP] Setting up distributions for evolution.");
			setupCoefficients();
			std::println("[DGLAP] Finished setting up distributions for evolution.");

			// if the next mass is zero, we are already done
			if (_alpha_s.masses(_nf+1) == 0.0)
			{
				std::println("[DGLAP] Next mass is zero. Quitting...");
				break;
			}

			
			// update all values
			_alpha_s.update(_nf);
			SplittingFunction::update(_nf, _alpha_s.beta0());
			_alpha0 = _alpha_s.post(_nf);
			_alpha1 = _alpha_s.pre(_nf+1);
			

			// determine if we are evolving to a
			// tabulated energy or threshold energy
			bool resum_tab = _qct<Qtab.size() && Qtab.at(_qct)<=_alpha_s.masses(_nf+1);
			bool resum_threshold = !resum_tab;
			
			double Q = Qtab[_qct]; // can do this no matter what
			
			// alpha1 needs to be manually calculated
			// if we are evolving to a tabulated energy
			if (resum_tab)
			{
				_alpha1 = _alpha_s.evaluate(_alpha_s.masses(_nf), Q, _alpha0);
				arr = std::ref(_F2);
				arr_singlet = std::ref(_F2);
			}
			else
			{
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
			if (_order == 1)
			{
				L2 = std::log((_alpha1*beta1 + 4.0*PI*beta0)
								/(_alpha0*beta1 + 4.0*PI*beta0));
			}
			else if (_order == 2)
			{
				L2 = std::log((16.0*PI_2*beta0 + 4.*PI*_alpha1*beta1 + _alpha1*_alpha1*beta2)
							/(16.*PI_2*beta0 + 4.*PI*_alpha0*beta1 + _alpha0*_alpha0*beta2));
				
				// analytic continuation for arctan
				double aux=4.0*beta0*beta2 - beta1*beta1;
				if (aux>=0)
				{
					L3 = std::atan(
						2.0*PI*(_alpha1-_alpha0)*std::sqrt(aux)
						/(2.*PI*(8.*PI*beta0+(_alpha1+_alpha0))+_alpha1*_alpha0*beta2)
					)/std::sqrt(aux);
				}
				else
				{
					L3 = std::tanh(
						2.*PI*(_alpha1-_alpha0)*std::sqrt(-aux)
						/(2.*PI*(8.*PI*beta0+(_alpha1+_alpha0))+_alpha1*_alpha0*beta2)
					)/std::sqrt(-aux);
				}
			}
			else if (_order == 3)
			{
				// NOTE: equation 2 and recursion relation 2 are determined from the log with the quadratic terms,
				// which were initially called L3, so I swapped them
				L3 = std::log((_alpha1 - r1)/(_alpha0 - r1));
				L2 = std::log((_alpha1*_alpha1 + b*_alpha1 + c)/(_alpha0*_alpha0 + b*_alpha0 + c));
				double aux = std::sqrt(-b*b + 4*c); // never negative, no need for analytic continuation
				L4 = std::atan((_alpha1-_alpha0)*aux / (2.0*_alpha0*_alpha1 + (_alpha0+_alpha1)*b + 2.0*c))/aux;
			}
			

			std::println("[DGLAP] Doing {} resummation", (resum_tab ? "tabulated" : "threshold" ));
			std::println("[DGLAP] AlphaS: {} --> {}", _alpha0, _alpha1);

			// only do evolution if alphas are different
			// (i.e. energy scales are different)
			if (_alpha0 != _alpha1)
			{
				std::println("[DGLAP] Starting singlet evolution and resummation...");
				evolveSinglet(arr_singlet, L1);
				std::println("[DGLAP] Finished singlet evolution and resummation.");

				std::println("[DGLAP] Starting non-singlet evolution and resummation...");
				_multi_thread ? 
					evolveNonSingletThreaded(arr, L1, L2, L3, L4) : 
					evolveNonSinglet(arr, L1, L2, L3, L4);
				std::println("[DGLAP] Finished non-singlet evolution and resummation.");

				std::println("[DGLAP] Fixing distributions...");
				fixDistributions(resum_tab, resum_threshold, temp_arr, temp_arr_singlet);
				std::println("[DGLAP] Finished fixing distributions.");

				// if we just resummed to a tabulated value,
				// _F2 contains our final distributions
				// we can just copy
				if (resum_tab)
				{
					final_dists = _F2;
				}
				else if (resum_threshold)
				{
					// if we just resummed to a threshold energy,
					// then we need to recopy the resultant distributions
					// from the temporary array
					// back to the n=0 piece
					for (uint j=0; j<DISTS; ++j)
					{
						switch (_order)
						{
							case 0: _A2[j][0] 		   = temp_arr[j]; break;
							case 1: _B2[j][0][0] 	   = temp_arr[j]; break;
							case 2: _C2[j][0][0][0]    = temp_arr[j]; break;
							case 3: _D2[j][0][0][0][0] = temp_arr[j]; break;
						}
						
					}
					for (uint j=0; j<=1; ++j)
						_S2[0][j][0] = temp_arr_singlet[j*31];
				}
			} // if (alpha0 != alpha1)

			if (resum_threshold && _order>=2 && _alpha_s.masses(_nf+2)!=0.0)
				heavyFlavorTreatment();
			
		} // for (_nf=_nfi; ; _nf++)

		std::println("[DGLAP] Done!");
		return final_dists;
	} // Evolve()

	
} // namespace Candia2
