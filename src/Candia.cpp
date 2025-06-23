#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/SplittingFn.hpp"

#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>


namespace Candia2
{

	DGLAPSolver::DGLAPSolver(const uint order, Grid const& grid, const double Qf,
							 std::shared_ptr<Distribution> initial_dist)
		: _order(order),  _grid(grid), _Qf(Qf),
		  _alpha_s(order, initial_dist->Q0(), initial_dist->Alpha0(), initial_dist->Masses()),
		  _dist(initial_dist)
	{	
		// reserve space in all the coefficient arrays
		std::cerr << "[DGLAP] Reserving space in coefficient arrays...\n";

		switch(_order)
		{
			case 0:
			{
				_A.resize(DISTS);
				for (uint j=0; j<DISTS; j++)
				{
					_A[j].resize(ITERATIONS);
					for (uint n=0; n<ITERATIONS; n++)
						_A[j][n].resize(_grid.Size());
				}
			}; break;
			case 1:
			{
				_B.resize(DISTS);
				for (uint j=0; j<DISTS; j++)
				{
					_B[j].resize(ITERATIONS);
					for (uint n=0; n<ITERATIONS; n++)
					{
						_B[j][n].resize(ITERATIONS);
						for (uint m=0; m<ITERATIONS; m++)
							_B[j][n][m].resize(_grid.Size());
					}
				}
			} break;
			case 2:
			{
				_C.resize(DISTS);
				for (uint j=0; j<DISTS; j++)
				{
					_C[j].resize(ITERATIONS);
					for (uint n=0; n<ITERATIONS; n++)
					{
						_C[j][n].resize(ITERATIONS);
						for (uint m=0; m<ITERATIONS; m++)
						{
							_C[j][n][m].resize(ITERATIONS);
							for (uint t=0; t<ITERATIONS; t++)
								_C[j][n][m][t].resize(_grid.Size());
						}
					}
				}
			} break;
			case 3:
			{
				_D.resize(DISTS);
				for (uint j=0; j<DISTS; j++)
				{
					_D[j].resize(ITERATIONS);
					for (uint n=0; n<ITERATIONS; n++)
					{
						_D[j][n].resize(ITERATIONS);
						for (uint m=0; m<ITERATIONS; m++)
						{
							_D[j][n][m].resize(ITERATIONS);
							for (uint t=0; t<ITERATIONS; t++)
							{
								_D[j][n][m][t].resize(ITERATIONS);
								for (uint s=0; s<ITERATIONS; s++)
									_D[j][n][m][t][s].resize(_grid.Size());
							}
						}
					}
				}
			} break;
			default: throw("Order " + std::to_string(_order) + " is invalid.");
		}

		std::cerr << "[DGLAP] Finished reserving non-singlet coefficients.\n";
		
		// singlet
		_S.resize(TRUNC_IDX+1);
		for (uint t=0; t<(TRUNC_IDX+1); t++)
		{
			_S[t].resize(2);
			for (uint j=0; j<2; j++)
			{
				_S[t][j].resize(ITERATIONS);
				for (uint n=0; n<ITERATIONS; n++)
					_S[t][j][n].resize(_grid.Size());
			}
		}

		std::cerr << "[DGLAP] Finished reserving singlet coefficients.\n";

		// final distributions
		_F.resize(DISTS);
		for (uint j=0; j<DISTS; j++)
			_F[j].resize(_grid.Size());

		
		std::cerr << "[DGLAP] Done allocating for all coefficients and the final distribution.\n";
		
		_alpha_s.CalculateThresholdValues(Qf);
		

		std::cerr << "[DGLAP] Creating splitting function pointers...\n";

		// just grab all of them, stores nearly nothing,
		// so it can't hurt
		_P0ns = std::make_shared<P0ns>();
		_P0qq = std::make_shared<P0qq>();
		_P0qg = std::make_shared<P0qg>();
		_P0gq = std::make_shared<P0gq>();	
		_P0gg = std::make_shared<P0gg>(_alpha_s.Beta0());

		_P1nsp = std::make_shared<P1nsp>();
		_P1nsm = std::make_shared<P1nsm>();
		_P1qq  = std::make_shared<P1qq>();
		_P1qg  = std::make_shared<P1qg>();
		_P1gq  = std::make_shared<P1gq>();	
		_P1gg  = std::make_shared<P1gg>();
			
		_P2nsp = std::make_shared<P2nsp>();
		_P2nsm = std::make_shared<P2nsm>();
		_P2nsv = std::make_shared<P2nsv>();
		_P2qq  = std::make_shared<P2qq>();
		_P2qg  = std::make_shared<P2qg>();
		_P2gq  = std::make_shared<P2gq>();	
		_P2gg  = std::make_shared<P2gg>();
			
		_P3nsp = std::make_shared<P3nsp>();
		_P3nsm = std::make_shared<P3nsm>();
		_P3nss = std::make_shared<P3nss>();

		std::cerr << "[DGLAP] Done creating splitting function pointers.\n";

		SetInitialConditions();

	    _nfi = _dist->Nfi();
		_nff = 6;
		while (_dist->Masses(_nff) > _Qf)
		{
			_nff--;
		}
	}

	DGLAPSolver::~DGLAPSolver()
	{
		std::cerr << "[DGLAP] Exiting... \n";
	}



	void DGLAPSolver::SetInitialConditions()
	{
		std::cerr << "[DGLAP] Setting initial conditions...\n";
	   
		// singlet
		for (uint k=0; k<_grid.Size(); k++)
		{
			double x = _grid[k];
			_S[0][0][0][k] = _dist->xg(x);
			_S[0][1][0][k] = _dist->xuv(x) + 2.0*_dist->xub(x)
				+ _dist->xdv(x) + 2.0*_dist->xdb(x) + 2.0*_dist->xs(x);
		}
	    
		switch (_order)
		{
			case 0:
			{
				for (uint k=0; k<_grid.Size(); k++)
				{
					double x = _grid[k];
					_A[7][0][k] = _dist->xub(x);
					_A[1][0][k] = _dist->xuv(x) + _A[7][0][k];
					_A[8][0][k] = _dist->xdb(x);
					_A[2][0][k] = _dist->xdv(x) + _A[8][0][k];
					_A[3][0][k] = _A[9][0][k] = _dist->xs(x);
				}
			} break;
			case 1:
			{
				for (uint k=0; k<_grid.Size(); k++)
				{
					double x = _grid[k];
					_B[7][0][0][k] = _dist->xub(x);
					_B[1][0][0][k] = _dist->xuv(x) + _B[7][0][0][k];
					_B[8][0][0][k] = _dist->xdb(x);
					_B[2][0][0][k] = _dist->xdv(x) + _B[8][0][0][k];
					_B[3][0][0][k] = _B[9][0][0][k] = _dist->xs(x);
				}
			} break;
			case 2:
			{
				for (uint k=0; k<_grid.Size(); k++)
				{
					double x = _grid[k];
					_C[7][0][0][0][k] = _dist->xub(x);
					_C[1][0][0][0][k] = _dist->xuv(x) + _C[7][0][0][0][k];
					_C[8][0][0][0][k] = _dist->xdb(x);
					_C[2][0][0][0][k] = _dist->xdv(x) + _C[8][0][0][0][k];
					_C[3][0][0][0][k] = _C[9][0][0][0][k] = _dist->xs(x);
				}
			} break;
			case 3:
			{
				for (uint k=0; k<_grid.Size(); k++)
				{
					double x = _grid[k];
					_D[7][0][0][0][0][k] = _dist->xub(x);
					_D[1][0][0][0][0][k] = _dist->xuv(x) + _D[7][0][0][0][0][k];
					_D[8][0][0][0][0][k] = _dist->xdb(x);
					_D[2][0][0][0][0][k] = _dist->xdv(x) + _D[8][0][0][0][0][k];
					_D[3][0][0][0][0][k] = _D[9][0][0][0][0][k] = _dist->xs(x);
				}

			} break;
		}

		std::cerr << "[DGLAP] Done setting initial conditions.\n";
	}



	double DGLAPSolver::RecRelS_1(std::vector<double> const& S, uint k,
								  std::shared_ptr<SplittingFunction> P)
	{
		return -(2.0)/_alpha_s.Beta0() * _grid.Convolution(S, P, k);
	}

	double DGLAPSolver::RecRelS_2(std::vector<double> const& S1,
								  std::vector<double> const& S2,
								  uint k,
								  std::shared_ptr<SplittingFunction> P0,
								  std::shared_ptr<SplittingFunction> P1)
	{
		
		double fact1 = -(2.0)/_alpha_s.Beta0() * _grid.Convolution(S1, P0, k);
		return (-fact1 - (1.0/(M_PI*_alpha_s.Beta0()))*_grid.Convolution(S2, P1, k));
	}



	double DGLAPSolver::RecRelLO(std::vector<double> const& A, uint k,
					  std::shared_ptr<SplittingFunction> P)
	{
		return -(2.0)/_alpha_s.Beta0() * _grid.Convolution(A, P, k);
	}

	double DGLAPSolver::RecRelNLO_1(std::vector<double> const& B,
									uint k,
									std::shared_ptr<SplittingFunction> P)
	{
		return -(2.0)/_alpha_s.Beta0() * _grid.Convolution(B, P, k);
	}

	double DGLAPSolver::RecRelNLO_2(std::vector<double> const& B1,
									std::vector<double> const& B2,
									uint k,
									std::shared_ptr<SplittingFunction> P)
	{
		return -B1[k] - (4.0/_alpha_s.Beta1() * _grid.Convolution(B2, P, k));
	}

	

	void DGLAPSolver::Evolve()
	{
		for (_nf=_nfi; ; _nf++)
		{
			
			std::cerr << "[DGLAP] Setting nf=" << _nf << '\n';
			SetupCoefficients();

			// if the next mass is zero, we are already done
			if (_alpha_s.Masses(_nf+1) == 0.0)
				break;

			// update all values
			_alpha_s.Update(_nf);
			SplittingFunction::UpdateNf(_nf);
			_alpha0 = _alpha_s.Post(_nf);
			_alpha1 = _alpha_s.Pre(_nf+1);
			std::cerr << "[DGLAP] Alpha_s: " << _alpha0 << " --> " << _alpha1 << '\n';

			_qct = 0;
			
			// only do evolution if alphas are different (i.e. energy scales are different
			if (_alpha0 != _alpha1)
			{
				// singlet sector
				EvolveSinglet();
				std::cerr << "[DGLAP] Evolve(): finished singlet evolution\n";
						
				// non-singlet sector
				EvolveNonSinglet();
				std::cerr << "[DGLAP] Evolve(): finished non-singlet evolution\n";

				// perform the resummation to the tabulated energy value
				Resum();
				std::cerr << "[DGLAP] Evolve(): finished resummation evolution\n";

				// finish resumming to the energy threshold
				ResumThreshold();
				std::cerr << "[DGLAP] Evolve(): finished resummation to threshold\n";
			}
		}
		// resetup the main quark/gluon distributions
		SetupFinalDistributions();

		std::cerr << "[DGLAP] Evolve(): Done!\n";
	}



	void DGLAPSolver::OutputDataFileNew(std::string const& filepath)
	{
		std::ofstream file(filepath);
		if (!file)
		{
			std::cerr << "[DGLAP] OutputDataFile(): Failure opening file '" << filepath << "'\n";
			exit(1);
		}
		
		file << std::fixed;
		file << std::setw(15) << std::setprecision(9);
		for (uint k=0; k<_grid.Size(); k++)
		{
			file << _grid.At(k);
			for (uint j=0; j<=12; j++)
				file << "  " << _F[j][k];
			file << '\n';
		}

		file.close();

		std::cerr << "[DGLAP] OutputDataFile(): Successfully saved data to file '" << filepath << "'\n";
	}

	void DGLAPSolver::OutputDataFileOld(std::string const& filepath)
	{
		FILE* outfile = fopen(filepath.c_str(), "w");
		for (uint k=0; k<_grid.Size(); k++) {
			fprintf(outfile, "%15.9lf", _grid.At(k));
			for (uint j=0; j<=12; j++)
				fprintf(outfile, "  %15.9lf", _F[j][k]);
			fprintf(outfile, "\n");
		}
		fclose(outfile);
	}


	

	void DGLAPSolver::SetupCoefficients()
	{
		switch (_order)
		{
			case 0: // LO
			{
				for (uint k=0; k<_grid.Size(); k++)
				{
					for (uint j=13; j<=18; j++)
						_A[j][0][k] = _A[j-12][0][k]-_A[j-6][0][k];
			
					_A[25][0][k]=0.;
					for (uint j=13; j<=18; j++)
						_A[25][0][k] += _A[j][0][k];
					for (uint j=26; j<=30; j++)
						_A[j][0][k] = _A[13][0][k]-_A[j-12][0][k];
					for (uint j=19; j<=24; j++)
						_A[j][0][k] = _A[j-18][0][k]+_A[j-12][0][k];

					_S[0][1][0][k] = 0.0;
					for (uint j=19; j<=24; j++)
						_S[0][1][0][k] += _A[j][0][k];
			
					for (uint j=32; j<=36; j++)
						_A[j][0][k]=_A[19][0][k]-_A[j-12][0][k];
				}
			} break;
			case 1: // NLO
			{
				for (uint k=0; k<_grid.Size(); k++)
				{
					for (uint j=13; j<=18; j++)
						_B[j][0][0][k] = _B[j-12][0][0][k]-_B[j-6][0][0][k];
			
					_B[25][0][0][k]=0.;
					for (uint j=13; j<=18; j++)
						_B[25][0][0][k] += _B[j][0][0][k];
					for (uint j=26; j<=30; j++)
						_B[j][0][0][k] = _B[13][0][0][k]-_B[j-12][0][0][k];
					for (uint j=19; j<=24; j++)
						_B[j][0][0][k] = _B[j-18][0][0][k]+_B[j-12][0][0][k];

					_S[0][1][0][k] = 0.0;
					for (uint j=19; j<=24; j++)
						_S[0][1][0][k] += _B[j][0][0][k];
			
					for (uint j=32; j<=36; j++)
						_B[j][0][0][k] = _B[19][0][0][k]-_B[j-12][0][0][k];
				}
			} break;
			case 2: // NNLO
			{
				for (uint k=0; k<_grid.Size(); k++)
				{
					for (uint j=13; j<=18; j++)
						_C[j][0][0][0][k] = _C[j-12][0][0][0][k]-_C[j-6][0][0][0][k];
			
					_C[25][0][0][0][k]=0.;
					for (uint j=13; j<=18; j++)
						_C[25][0][0][0][k] += _C[j][0][0][0][k];
					for (uint j=26; j<=30; j++)
						_C[j][0][0][0][k] = _C[13][0][0][0][k]-_C[j-12][0][0][0][k];
					for (uint j=19; j<=24; j++)
						_C[j][0][0][0][k] = _C[j-18][0][0][0][k]+_C[j-12][0][0][0][k];

					_S[0][1][0][k] = 0.0;
					for (uint j=19; j<=24; j++)
						_S[0][1][0][k] += _C[j][0][0][0][k];
			
					for (uint j=32; j<=36; j++)
						_C[j][0][0][0][k] = _C[19][0][0][0][k]-_C[j-12][0][0][0][k];
				}
			} break;
			case 3: //N3LO
			{
				for (uint k=0; k<_grid.Size(); k++)
				{
					for (uint j=13; j<=18; j++)
						_D[j][0][0][0][0][k] = _D[j-12][0][0][0][0][k]-_D[j-6][0][0][0][0][k];
			
					_D[25][0][0][0][0][k]=0.;
					for (uint j=13; j<=18; j++)
						_D[25][0][0][0][0][k] += _D[j][0][0][0][0][k];
					for (uint j=26; j<=30; j++)
						_D[j][0][0][0][0][k] = _D[13][0][0][0][0][k]-_D[j-12][0][0][0][0][k];
					for (uint j=19; j<=24; j++)
						_D[j][0][0][0][0][k] = _D[j-18][0][0][0][0][k]+_D[j-12][0][0][0][0][k];

					_S[0][1][0][k] = 0.0;
					for (uint j=19; j<=24; j++)
						_S[0][1][0][k] += _D[j][0][0][0][0][k];
			
					for (uint j=32; j<=36; j++)
						_D[j][0][0][0][0][k] = _D[19][0][0][0][0][k]-_D[j-12][0][0][0][0][k];
				}
			} break;
			default: throw("Order " + std::to_string(_order) + " is invalid");
		}
	}

	void DGLAPSolver::SetupFinalDistributions()
	{
	    // reconstruct the actual quark distributions
		// for tabulated Q values
		double Nf = static_cast<double>(_nf);
		
		for (uint k=0; k<_grid.Size()-1; k++)
		{
			_F[19][k] = _F[31][k];

			_F[19][k] = _F[31][k];
			for (uint j=32; j<=30+_nf; j++)
				_F[19][k] += _F[j][k];
			_F[19][k] /= Nf;

			for (uint j=20; j<=18+_nf; j++)
				_F[j][k] = _F[19][k] - _F[j+12][k];

			for (uint j=1; j<=_nf; j++) {
				_F[j][k]   = 0.5*(_F[j+18][k] + _F[j+12][k]);
				_F[j+6][k] = 0.5*(_F[j+18][k] - _F[j+12][k]);
			}

			_F[25][k]=0.0;
			for (uint j=13; j<=12+_nf; j++)
				_F[25][k] += _F[j][k];
			for (uint j=26; j<=24+_nf; j++)
				_F[j][k] = _F[13][k] - _F[j-12][k];
		}
	}



	void DGLAPSolver::EvolveSinglet()
	{
		for (uint n=0; n<ITERATIONS-1; n++)
		{
			std::cerr << "[DGLAP] EvolveSinglet(): Iteration " << n << '\n';

			// singlet 0-coefficients
			for (uint k=0; k<_grid.Size()-1; k++)
			{
				_S[0][1][n+1][k] = RecRelS_1(_S[0][1][n], k, _P0qq) + RecRelS_1(_S[0][0][n], k, _P0qg);
				_S[0][0][n+1][k] = RecRelS_1(_S[0][1][n], k, _P0gq) + RecRelS_1(_S[0][0][n], k, _P0gg);
			}

			if (_order > 0)
			{
				// offset of singlet 1-coefficients
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					for (uint j=0; j<2; j++)
						_S[1][j][n+1][k] = -_S[0][j][n+1][k] * _alpha_s.Beta1()/(4.0*M_PI*_alpha_s.Beta0()) - _S[1][j][n][k];
				}

				// evolution of singlet 1-coefficients
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					_S[1][1][n+1][k] += RecRelS_2(_S[0][1][n], _S[1][1][n], k, _P0qq, _P1qq) + RecRelS_2(_S[0][0][n], _S[1][0][n], k, _P0qg, _P1qg);
					_S[1][0][n+1][k] += RecRelS_2(_S[0][1][n], _S[1][1][n], k, _P0gq, _P1gq) + RecRelS_2(_S[0][0][n], _S[1][0][n], k, _P0gg, _P1gg);
				}
			}

			// further evolution based on truncation index
			for (uint t=2; t<TRUNC_IDX+1; t++)
			{
				double T = static_cast<double>(t);
				// singlet n-coefficients (after doing 0 and 1)
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					for (uint j=0; j<2; j++)
					{
						_S[t][j][n+1][k] = -_S[t-1][j][n+1][k] * _alpha_s.Beta1()/(4.0*M_PI*_alpha_s.Beta0())
							- T*_S[t][j][n][k] - (T-1.0)*_S[t-1][j][n][k] * _alpha_s.Beta1()/(4.0*M_PI*_alpha_s.Beta0());

						if (_order == 2)
						{
							_S[t][j][n+1][k] -= _S[t-2][j][n+1][k] * _alpha_s.Beta2()/(16.0*M_PI_2*_alpha_s.Beta0())
								+ (T-2.0)*_S[t-2][j][n][k] * _alpha_s.Beta2()/(16.0*M_PI_2*_alpha_s.Beta0());
						}
					}
				}

				// NLO singlet n-coefficients
				if (_order == 1)
				{
					for (uint k=0; k<_grid.Size()-1; k++)
					{
						_S[t][1][n+1][k] += RecRelS_2(_S[t-1][1][n], _S[t][1][n], k, _P0qq, _P1qq)
							+ RecRelS_2(_S[t-1][0][n], _S[t][0][n], k, _P0qg, _P1qg);
						_S[t][0][n+1][k] += RecRelS_2(_S[t-1][1][n], _S[t][1][n], k, _P0gq, _P1gq)
							+ RecRelS_2(_S[t-1][0][n], _S[t][0][n], k, _P0gg, _P1gg);
					}
				}
			}
		}
	}




	void DGLAPSolver::EvolveNonSinglet()
	{
		switch (_order)
		{
			case 0: // LO
			{
				// solve for all A_n(x) coefficients
				for (uint n=0; n<ITERATIONS-1; n++)
				{
					std::cerr << "[DGLAP] EvolveNonSinglet() LO Iteration " << n << '\n';

					for (uint k=0; k<_grid.Size()-1; k++)
					{
						for (uint j=13; j<=12; j++)
							_A[j][n+1][k] = RecRelLO(_A[j][n], k, _P0ns);
					
						for (uint j=32; j<=30+_nf; j++)
							_A[j][n+1][k] = RecRelLO(_A[j][n], k, _P0ns);
					}
				}
			} break;
			case 1: // NLO
			{
				for (uint s=1; s<ITERATIONS; s++)
				{
					std::cerr << "[DGLAP] EvolveNonSinglet() NLO Iteration " << s << '\n';

					for (uint k=0; k<_grid.Size()-1;k++)
					{
						for (uint j=13; j<=12+_nf; j++)
						{
							for (uint n=1; n<=s; n++)
								_B[j][s][n][k] = RecRelNLO_1(_B[j][s-1][n-1], k, _P0ns);

							_B[j][s][0][k] = RecRelNLO_2(_B[j][s][0], _B[j][s-1][0], k, _P1nsm);
						}

						for (uint j=32; j<=30+_nf; j++)
						{
							for (uint n=1; n<=s; n++)
								_B[j][s][n][k] = RecRelNLO_1(_B[j][s-1][n-1], k, _P0ns);

							_B[j][s][0][k] = RecRelNLO_2(_B[j][s][1], _B[j][s-1][0], k, _P1nsp);
						}
					}
				}
				
			} break;
			case 2: // NNLO
			{
				std::cerr << "[DGLAP] EvolveNonSinglet() NNLO Evolution not implemented yet.\n";
				exit(1);
			} break;
			case 3: // N3LO
			{
				std::cerr << "[DGLAP] EvolveNonSinglet() N3LO Evolution not implemented yet.\n";
				exit(1);
			} break;
			default:
			{
				std::cerr << "[DGLAP] EvolveNonSinglet() Encountered invalid order: " << _order << '\n';
				exit(1);
			}
		}
	}



	void DGLAPSolver::Resum()
	{
		std::cerr << "[DGLAP] Resum(): resumming to tabulated energy\n";
		
		std::vector<double> Qtab{ _Qf };
		
		for ( ; Qtab[_qct]<=_alpha_s.Masses(_nf+1) && _qct<Qtab.size(); _qct++)
		{
			_alpha1 = _alpha_s.Evaluate(_alpha_s.Masses(_nf), Qtab[_qct], _alpha0);
			double L1 = std::log(_alpha1/_alpha0);
			double L2  = std::log((_alpha1*_alpha_s.Beta1() + 4.0*M_PI*_alpha_s.Beta0())
								  /(_alpha0*_alpha_s.Beta1() + 4.0*M_PI*_alpha_s.Beta0()));


			std::cerr << "[DGLAP] Resum(): _alpha1 = " << _alpha1 << '\n';
			std::cerr << "[DGLAP] Resum(): L1=" << L1 << ", L2=" << L2 << '\n';

			// singlet
			for (uint j=0; j<=1; j++)
			{
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					_F[j*31][k] = _S[0][j][0][k];
						
					for (uint n=0; n<ITERATIONS; n++)
					{
						for (uint t=0; t<TRUNC_IDX+1; t++)
							_F[j*31][k] += _S[t][j][n][k] * std::pow(_alpha1, t)*std::pow(L1, n)/Factorial(n);
						
					}
				}
			}


			// non-singlet
			switch (_order)
			{
				case 0: // LO
				{
					for (uint k=0; k<_grid.Size()-1; k++)
					{
						for (uint j=13; j<=30+_nf; j++)
						{
							_F[j][k] = _A[j][0][k];

							if (j == (12+_nf))
								j = 31;
						}
						
						for (uint n=1; n<ITERATIONS; n++)
						{
							for (uint j=13; j<=12+_nf; j++)
								_F[j][k] += _A[j][n][k]*std::pow(L1, n)/Factorial(n);
							for (uint j=32; j<=30+_nf; j++)
								_F[j][k] += _A[j][n][k]*std::pow(L1, n)/Factorial(n);
						}
					}
				}; break;
				case 1: // NLO
				{
					for (uint k=0; k<_grid.Size()-1; k++)
					{
						for (uint j=13; j<=30+_nf; j++)
						{
							_F[j][k] = _B[j][0][0][k];

							if (j==12+_nf)
								j=31;
						}

						for (uint s=1; s<ITERATIONS; s++)
						{
							for (uint n=0; n<=s; n++)
							{
								for (uint j=13; j<=12+_nf; j++)
									_F[j][k]+=_B[j][s][n][k]/Factorial(n)/Factorial(s-n)*std::pow(L1,n)*std::pow(L2,(s-n));

								for (uint j=32; j<=30+_nf; j++)
									_F[j][k]+=_B[j][s][n][k]/Factorial(n)/Factorial(s-n)*std::pow(L1,n)*std::pow(L2,(s-n));
							}
						}

						
					}
				}; break;
			}
		}
	}



	void DGLAPSolver::ResumThreshold()
	{
		// now sum the recursive series for threshold Q values i guess?
		_alpha1 = _alpha_s.Pre(_nf+1);
		double L1 = std::log(_alpha1/_alpha0);
		double L2 = std::log((_alpha1*_alpha_s.Beta1() + 4.0*M_PI*_alpha_s.Beta0())
							 /(_alpha0*_alpha_s.Beta1() + 4.0*M_PI*_alpha_s.Beta0()));
		double Nf = static_cast<double>(_nf);


		// singlet 
		for (uint j=0; j<2; j++)
		{
			for (uint k=0; k<_grid.Size()-1; k++)
			{
				for (uint n=1; n<ITERATIONS; n++)
				{
					for (uint t=0; t<TRUNC_IDX+1; t++)
						_S[0][j][0][k] += _S[t][j][n][k]*std::pow(_alpha1, t)*std::pow(L1, n)/Factorial(n);
				}
			}
		}


		// non-singlet
		switch (_order)
		{
			case 0: // LO
			{
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					// resum
					for (uint n=1; n<ITERATIONS; n++)
					{
						for (uint j=12; j<=12+_nf; j++)
							_A[j][0][k] += _A[j][n][k]*pow(L1, n)/Factorial(n);
						for (uint j=32; j<=30+_nf; j++)
							_A[j][0][k] += _A[j][n][k]*pow(L1, n)/Factorial(n);
					}
				}

				// set distributions
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					_A[19][0][k] = _S[0][1][0][k];
					for (uint j=32; j<=30+_nf; j++)
						_A[19][0][k] += _A[j][0][k];
					_A[19][0][k] /= Nf;

					for (uint j=20; j<=18+_nf; j++)
						_A[j][0][k] = _A[19][0][k] - _A[j+12][0][k];

					for (uint j=1; j<=_nf; j++)
					{
						_A[j][0][k]   = 0.5*(_A[j+18][0][k] + _A[j+12][0][k]);
						_A[j+6][0][k] = 0.5*(_A[j+18][0][k] - _A[j+12][0][k]);
					}
				}
			} break;
			case 1: // NLO
			{
				// resum
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					for (uint s=1; s<ITERATIONS; s++)
					{
						for (uint n=0; n<=s; n++)
						{
							for (uint j=13; j<=12+_nf; j++)
								_B[j][0][0][k] += _B[j][s][n][k]/Factorial(n)/Factorial(s-n)*std::pow(L1,n)*std::pow(L2,(s-n));

							for (uint j=32; j<=30+_nf; j++)
								_B[j][0][0][k] += _B[j][s][n][k]/Factorial(n)/Factorial(s-n)*std::pow(L1,n)*std::pow(L2,(s-n));
						}
					}
				}

				// set distributions
				for (uint k=0; k<_grid.Size()-1;k++) {
					_B[19][0][0][k]=_S[0][1][0][k];
					for (uint j=32; j<=30+_nf; j++)
						_B[19][0][0][k] += _B[j][0][0][k];
					_B[19][0][0][k]/=Nf;

					for (uint j=20; j<=18+_nf; j++)
						_B[j][0][0][k] = _B[19][0][0][k]-_B[j+12][0][0][k];

					for (uint j=1; j<=_nf; j++)
					{
						_B[j][0][0][k]   = 0.5*(_B[j+18][0][0][k]+_B[j+12][0][0][k]);
						_B[j+6][0][0][k] = 0.5*(_B[j+18][0][0][k]-_B[j+12][0][0][k]);
					}
				}
			}; break;
		}
	}

	
	
}
