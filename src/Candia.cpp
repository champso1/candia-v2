#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/Math.hpp"
#include "Candia-v2/SplittingFn.hpp"


#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

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
			} break;
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

				// also initialize the N3LO coefficients we will need
				_r1[3] = -0.0884922; _b[3] = -2.0 * 0.0176202; _c[3] = std::pow(0.0176202, 2) + std::pow(0.090083, 2);
				_r1[4] = -0.0962106; _b[4] = -2.0 * 0.0228195; _c[4] = std::pow(0.0228195, 2) + std::pow(0.101286, 2);
				_r1[5] = -0.1050890; _b[5] = -2.0 * 0.0338021; _c[5] = std::pow(0.0338021, 2) + std::pow(0.118211, 2);
				_r1[6] = -0.1136210; _b[6] = -2.0 * 0.0633832; _c[6] = std::pow(0.0633832, 2) + std::pow(0.144577, 2);

				// technially the above values are the solutions
				// for a_s, NOT alpha_s
				// this therefore makes r1, r2, and r3 a factor of 4pi off
				// b is 2*Re[r2], so it gets a single power of 4pi,
				// but c is |r2|^2 (or |r3|^2) so it gets (4pi)^2
				std::for_each(_r1.begin(), _r1.end(), [](double &x){
					x *= 4*PI;
				});
				std::for_each(_b.begin(), _b.end(), [](double &x){
					x *= 4*PI;
				});
				std::for_each(_c.begin(), _c.end(), [](double &x){
					x *= std::pow(4*PI, 2);
				});
			} break;
			default:
			{
				std::cerr << "[DGLAP] The order " << _order << " is invalid. Must be between 0 and 3 (inclusive).\n";
				exit(1);
			}
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

		
		std::cerr << "[DGLAP] Done allocating for all coefficients and the final distributions.\n";

		
		_alpha_s.CalculateThresholdValues(Qf);
		

		std::cerr << "[DGLAP] Creating splitting function pointers...\n";

		// just grab all of them, stores nearly nothing,
		// so it can't hurt
		_P0ns = std::make_shared<P0ns>();
		_P0qq = std::make_shared<P0qq>();
		_P0qg = std::make_shared<P0qg>();
		_P0gq = std::make_shared<P0gq>();	
		_P0gg = std::make_shared<P0gg>();

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
		_P3nsv = std::make_shared<P3nsv>();

		std::cerr << "[DGLAP] Done creating splitting function pointers.\n";

		FillSplittingFunctionCaches();

		std::cerr << "[DGLAP] Creating operator matrix element pointers...\n";

		_A2ns = std::make_shared<A2ns>();
		_A2gq = std::make_shared<A2gq>();
		_A2gg = std::make_shared<A2gg>();
		_A2hq = std::make_shared<A2hq>();
		_A2hg = std::make_shared<A2hg>();

		std::cerr << "[DGLAP] Done creating operator matrix element pointers.\n";

		SetInitialConditions();

	    _nfi = _dist->Nfi();
		_nff = _alpha_s.Nff(_nfi, _Qf);
		std::cerr << "[DGLAP] Evolving to " << _nff << " flavors.\n";

		_output_file_index = 1;
		_qct = 0;
	}

	DGLAPSolver::~DGLAPSolver()
	{
		std::cerr << "[DGLAP] Exiting... \n";
	}



	void DGLAPSolver::SetInitialConditions()
	{
		std::cerr << "[DGLAP] Setting initial conditions...\n";
	   
		// singlet
		for (uint k=0; k<_grid.Size()-1; k++)
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
				for (uint k=0; k<_grid.Size()-1; k++)
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
				for (uint k=0; k<_grid.Size()-1; k++)
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
				for (uint k=0; k<_grid.Size()-1; k++)
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
				for (uint k=0; k<_grid.Size()-1; k++)
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

	void DGLAPSolver::FillSplittingFunctionCaches()
	{
		std::cerr << "[DGLAP] FillSplittingFunctionCaches(): NOT YET IMPLEMENTED\n";
	}

	double DGLAPSolver::RecRelS_1(std::vector<double> const& S, uint k,
								  std::shared_ptr<SplittingFunction> P)
	{
		double conv = _grid.Convolution(S, P, k);
		return -(2.0/_alpha_s.Beta0()) * conv;
	}

	double DGLAPSolver::RecRelS_2(std::vector<double> const& S1,
								  std::vector<double> const& S2,
								  uint k,
								  std::shared_ptr<SplittingFunction> P0,
								  std::shared_ptr<SplittingFunction> P1)
	{
		double conv1 = _grid.Convolution(S2, P0, k);
		double conv2 = _grid.Convolution(S1, P1, k);
		double res = conv1*2.0/_alpha_s.Beta0();
		res += conv2 / (PI*_alpha_s.Beta0());
		
		return -res;
	}

	double DGLAPSolver::RecRelS_3(std::vector<double> const& S1,
								  std::vector<double> const& S2,
								  std::vector<double> const& S3,
								  uint k,
								  std::shared_ptr<SplittingFunction> P0,
								  std::shared_ptr<SplittingFunction> P1,
								  std::shared_ptr<SplittingFunction> P2)
	{
		double conv1 = _grid.Convolution(S3, P0, k);
		double conv2 = _grid.Convolution(S2, P1, k);
		double conv3 = _grid.Convolution(S1, P2, k);

		double res = conv1 * 2.0/_alpha_s.Beta0();
		res += conv2 / (PI*_alpha_s.Beta0());
		res += conv3 / (2.0*PI_2*_alpha_s.Beta0());
		
		return -res;
	}



	double DGLAPSolver::RecRelLO(std::vector<double> const& A, uint k,
					  std::shared_ptr<SplittingFunction> P)
	{
		double p1 = -2.0/_alpha_s.Beta0();
		double p2 = _grid.Convolution(A, P, k);
		return p1*p2;
	}


	double DGLAPSolver::RecRelNLO_1(std::vector<double> const& B,
									uint k,
									std::shared_ptr<SplittingFunction> P)
	{
		return -(2.0/_alpha_s.Beta0()) * _grid.Convolution(B, P, k);
	}

	double DGLAPSolver::RecRelNLO_2(std::vector<double> const& B,
									uint k,
									std::shared_ptr<SplittingFunction> P)
	{
		double res = -(4.0/_alpha_s.Beta1()) * _grid.Convolution(B, P, k);
		return res;
	}


	double DGLAPSolver::RecRelNNLO_1(std::vector<double> const& C,
									 uint k,
									 std::shared_ptr<SplittingFunction> P)
	{
		double conv = _grid.Convolution(C, P, k);
		return -(2.0)/_alpha_s.Beta0() * conv;
	}

	double DGLAPSolver::RecRelNNLO_2(std::vector<double> const& C,
									 uint k,
									 std::shared_ptr<SplittingFunction> P)
	{
		double conv = _grid.Convolution(C, P, k);
		double res = -4.0/_alpha_s.Beta2() * conv;
		return res;
	}

	double DGLAPSolver::RecRelNNLO_3(std::vector<double> const& C,
									 uint k,
									 std::shared_ptr<SplittingFunction> P)
	{
		double conv = _grid.Convolution(C, P, k);
		return -8.0 * conv;
	}


	double DGLAPSolver::RecRelN3LO_1(std::vector<double> const& D,
									 uint k,
									 std::shared_ptr<SplittingFunction> P)
	{
		double conv = _grid.Convolution(D, P, k);
		return -(2.0)/_alpha_s.Beta0() * conv;
	}

	double DGLAPSolver::RecRelN3LO_2(std::vector<double> const& D,
									 uint k,
									 std::shared_ptr<SplittingFunction> P1,
									 std::shared_ptr<SplittingFunction> P2,
									 std::shared_ptr<SplittingFunction> P3)
	{
		// for simplified notation
		const double r1 = _r1[_nf];
		// const double b = _b[_nf];
		// const double c = _c[_nf];

		double conv1 = _grid.Convolution(D, P1, k);
		double conv2 = _grid.Convolution(D, P2, k);
		double conv3 = _grid.Convolution(D, P3, k);

		const double fac1 = -4*conv1;
		const double fac2 = -8*r1*conv2;
		const double fac3 = -16*r1*r1*conv3;
		
	    return fac1+fac2+fac3;
	}

	
	double DGLAPSolver::RecRelN3LO_3(std::vector<double> const& D,
									 uint k,
									 std::shared_ptr<SplittingFunction> P1,
									 std::shared_ptr<SplittingFunction> P2,
									 std::shared_ptr<SplittingFunction> P3)
	{
		// for simplified notation
		const double r1 = _r1[_nf];
		const double b = _b[_nf];
		const double c = _c[_nf];

		double conv1 = _grid.Convolution(D, P1, k);
		double conv2 = _grid.Convolution(D, P2, k);
		double conv3 = _grid.Convolution(D, P3, k);

		const double fac1 = 2*conv1;
		const double fac2 = 4*r1*conv2;
		const double fac3 = -8*(c + b*r1)*conv3;
		
		double res = fac1+fac2+fac3;
	    return res;
	}

	double DGLAPSolver::RecRelN3LO_4(std::vector<double> const& D,
									 uint k,
									 std::shared_ptr<SplittingFunction> P1,
									 std::shared_ptr<SplittingFunction> P2,
									 std::shared_ptr<SplittingFunction> P3)
	{
		const double r1 = _r1[_nf];
		const double b = _b[_nf];
		const double c = _c[_nf];

		double conv1 = _grid.Convolution(D, P1, k);
		double conv2 = _grid.Convolution(D, P2, k);
		double conv3 = _grid.Convolution(D, P3, k);

		const double fac1 = 8*(b+r1)*conv1;
		const double fac2 = -16*c*conv2;
		const double fac3 = -32*c*r1*conv3;
		
	    return fac1+fac2+fac3;
	}


	

	void DGLAPSolver::Evolve()
	{
		std::vector<double> Qtab{ _Qf };

		for (_nf=_nfi; ; _nf++)
		{
			
			std::cerr << "[DGLAP] Setting nf=" << _nf << '\n';
			SetupCoefficients();

			// if the next mass is zero, we are already done
			if (_alpha_s.Masses(_nf+1) == 0.0)
				break;

			// update all values
			_alpha_s.Update(_nf);
			SplittingFunction::Update(_nf, _alpha_s.Beta0());
			_alpha0 = _alpha_s.Post(_nf);
			_alpha1 = _alpha_s.Pre(_nf+1);
			std::cerr << "[DGLAP] Alpha_s: " << _alpha0 << " --> " << _alpha1 << '\n';

			// only do evolution if alphas are different (i.e. energy scales are different
			if (_alpha0 != _alpha1)
			{
				// singlet sector
				EvolveSinglet();
				std::cerr << "[DGLAP] Evolve(): finished singlet evolution\n";
						
				// non-singlet sector
				EvolveNonSinglet();
				std::cerr << "[DGLAP] Evolve(): finished non-singlet evolution\n";

				// perform the resummation to the tabulated energy value(s)
				for ( ; Qtab[_qct]<=_alpha_s.Masses(_nf+1) && _qct<Qtab.size(); _qct++)
					Resum(Qtab[_qct]);
				std::cerr << "[DGLAP] Evolve(): finished resummation evolution\n";

				// finish resumming to the energy threshold
				ResumThreshold();
				std::cerr << "[DGLAP] Evolve(): finished resummation to threshold\n";
			}

			if (_order>=2 && _alpha_s.Masses(_nf+2)!=0.0)
				HeavyFlavorTreatment();
		}

		// resetup the main quark/gluon distributions
		// TODO: figure out why I had this here, as it is useless
		// and im pretty sure also just incorrect
		// SetupFinalDistributions();

		std::cerr << "[DGLAP] Evolve(): Done!\n";
	}



	void DGLAPSolver::SetOutputDataFileNew(std::string const& filepath)
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
	void DGLAPSolver::SetOutputDataFileOld(std::string const& filepath)
	{
		UNUSED(filepath);
	    return;
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
			
					_C[25][0][0][0][k]=0.0;
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
			}
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
			std::cerr << "[DGLAP] EvolveSinglet(): Singlet iteration " << n << '\n';

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
					for (uint j=0; j<=1; j++)
						_S[1][j][n+1][k] = -_S[0][j][n+1][k] * _alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0()) - _S[1][j][n][k];
				}

				// evolution of singlet 1-coefficients
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					_S[1][1][n+1][k] += RecRelS_2(_S[0][1][n], _S[1][1][n], k, _P0qq, _P1qq) + RecRelS_2(_S[0][0][n], _S[1][0][n], k, _P0qg, _P1qg);
					_S[1][0][n+1][k] += RecRelS_2(_S[0][1][n], _S[1][1][n], k, _P0gq, _P1gq) + RecRelS_2(_S[0][0][n], _S[1][0][n], k, _P0gg, _P1gg);
				}
			}

			// further evolution based on truncation index
			for (uint t=2; t<=TRUNC_IDX; t++)
			{
				double T = static_cast<double>(t);
				// singlet n-coefficients (after doing 0 and 1)
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					for (uint j=0; j<=1; j++)
					{
						_S[t][j][n+1][k] = -_S[t-1][j][n+1][k] * _alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0())
							- T*_S[t][j][n][k] - (T-1.0)*_S[t-1][j][n][k]*_alpha_s.Beta1()/(4.0*PI*_alpha_s.Beta0());

						if (_order >= 2)
						{
							_S[t][j][n+1][k] -= _S[t-2][j][n+1][k] * _alpha_s.Beta2()/(16.0*PI_2*_alpha_s.Beta0())
								+ (T-2.0)*_S[t-2][j][n][k] * _alpha_s.Beta2()/(16.0*PI_2*_alpha_s.Beta0());
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

				// NNLO (& N3LO) singlet coefficients
				if (_order >= 2)
				{
				    for (uint k=0; k<_grid.Size()-1; k++)
					{
						_S[t][1][n+1][k] += RecRelS_3(_S[t-2][1][n], _S[t-1][1][n], _S[t][1][n], k, _P0qq, _P1qq, _P2qq)
										  + RecRelS_3(_S[t-2][0][n], _S[t-1][0][n], _S[t][0][n], k, _P0qg, _P1qg, _P2qg);
						_S[t][0][n+1][k] += RecRelS_3(_S[t-2][1][n], _S[t-1][1][n], _S[t][1][n], k, _P0gq, _P1gq, _P2gq)
										  + RecRelS_3(_S[t-2][0][n], _S[t-1][0][n], _S[t][0][n], k, _P0gg, _P1gg, _P2gg);
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
					std::cerr << "[DGLAP] EvolveNonSinglet(): LO Iteration " << n << '\n';

					for (uint k=0; k<_grid.Size()-1; k++)
					{
						for (uint j=13; j<=12+_nf; j++)
							_A[j][n+1][k] = RecRelLO(_A[j][n], k, _P0ns);
					
						for (uint j=32; j<=30+_nf; j++)
							_A[j][n+1][k] = RecRelLO(_A[j][n], k, _P0ns);

					}
				}
			} break;
			case 1: // NLO
			{
				for (uint s=1; s<=ITERATIONS-1; s++)
				{
					std::cerr << "[DGLAP] EvolveNonSinglet() NLO Iteration " << s << '\n';

					for (uint k=0; k<_grid.Size()-1;k++)
					{
						for (uint j=13; j<=12+_nf; j++)
						{
							for (uint n=1; n<=s; n++)
								_B[j][s][n][k] = RecRelNLO_1(_B[j][s-1][n-1], k, _P0ns);

							_B[j][s][0][k] = -_B[j][s][1][k] + RecRelNLO_2(_B[j][s-1][0], k, _P1nsm);
						}

						for (uint j=32; j<=30+_nf; j++)
						{
							for (uint n=1; n<=s; n++)
								_B[j][s][n][k] = RecRelNLO_1(_B[j][s-1][n-1], k, _P0ns);

							_B[j][s][0][k] = -_B[j][s][1][k] + RecRelNLO_2(_B[j][s-1][0], k, _P1nsp);
						}
					}
				}
				
			} break;
			case 2: // NNLO
			{
				Clock::time_point t0, tf;

				for (uint s=1; s<=ITERATIONS-1; s++)
				{
					t0 = Clock::now();

					std::cerr << "[DGLAP] EvolveNonSinglet(): NNLO Iteration " << s << '\n';

					for (uint k=0; k<_grid.Size()-1; k++)
					{
						for (uint j=26; j<=24+_nf; j++)
						{
							for (uint t=1; t<=s; t++)
							{
								for (uint n=1; n<=t; n++)
									_C[j][s][t][n][k] = RecRelNNLO_1(_C[j][s-1][t-1][n-1], k, _P0ns);
							}

							double fac1 = -0.5*_C[j][s][s][1][k];
							double fac2 = RecRelNNLO_2(_C[j][s-1][s-1][0], k, _P2nsm);
							_C[j][s][s][0][k] = fac1 + fac2;

							// these must be regular ints;
							// unsigned ints, when they are 0 and get --,
							// underflow back to positive 4b,
							// remaining positive and the loop continues (very bad!)
							for (int t=s-1; t>=0; t--)
							{	
								double fac1 = -2.0*_alpha_s.Beta1()*(_C[j][s][t+1][0][k] + _C[j][s][t+1][1][k]);
								double fac2 = RecRelNNLO_3(_C[j][s-1][t][0], k, _P1nsm);
								_C[j][s][t][0][k] = fac1 + fac2;
							}
						}

						for (uint j=32; j<31+_nf; j++)
						{
							for (uint t=1;t<=s;t++)
							{
								for (uint n=1; n<=t; n++)
									_C[j][s][t][n][k] = RecRelNNLO_1(_C[j][s-1][t-1][n-1], k, _P0ns);
							}

							_C[j][s][s][0][k] = -0.5*_C[j][s][s][1][k] + RecRelNNLO_2(_C[j][s-1][s-1][0], k, _P2nsp);

							for (int t=s-1; t>=0; t--)
								_C[j][s][t][0][k] = -2.0*_alpha_s.Beta1()*(_C[j][s][t+1][0][k] + _C[j][s][t+1][1][k]) + RecRelNNLO_3(_C[j][s-1][t][0], k, _P1nsp);
						}

						for (uint t=1; t<=s; t++)
						{
							for (uint n=1; n<=t; n++)
								_C[25][s][t][n][k] = RecRelNNLO_1(_C[25][s-1][t-1][n-1], k, _P0ns);
						}

						_C[25][s][s][0][k] = -0.5*_C[25][s][s][1][k] + RecRelNNLO_2(_C[25][s-1][s-1][0], k, _P2nsv);

						for (int t=s-1; t>=0; t--)
							_C[25][s][t][0][k] = -2.0*_alpha_s.Beta1()*(_C[25][s][t+1][0][k] + _C[25][s][t+1][1][k]) + RecRelNNLO_3(_C[25][s-1][t][0], k, _P1nsm);
					}

					tf = Clock::now();
					std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(tf-t0);
					std::cerr << "[DGLAP] EvolveNonSinglet(): NNLO iteration took " << elapsed.count() << " seconds.\n";
				}
			} break;
			case 3: // N3LO
			{
				// some shorthand
				double r1 = _r1[_nf];
				double b = _b[_nf];
				double c = _c[_nf];
				double gamma = r1*r1 + r1*b + c;
				double fac = c + b*r1;

				Clock::time_point t0, tf;
				
				for (uint s=1; s<=ITERATIONS-1; s++)
				{
					t0 = Clock::now();
						
					std::cerr << "[DGLAP] EvolveNonSinglet(): N3LO iteration " << s << '\n';
					
					for (uint k=0; k<_grid.Size()-1; k++)
					{
						for (uint j=26; j<25+_nf; j++)
						{
							// RecRel #1:
							for (uint t=1; t<=s; t++)
							{
								for (uint m=1; m<=t; m++)
								{
									for (uint n=1; n<=m; n++)
										_D[j][s][t][m][n][k] = RecRelN3LO_1(_D[j][s-1][t-1][m-1][n-1], k, _P0ns);
								}
							}

							// RecRel #2:
							{
								double fac1 = 0.5*(_alpha_s.Beta1() + r1*_alpha_s.Beta2() - fac*_alpha_s.Beta3()) *_D[j][s][s][s][1][k];
								double fac2 = RecRelN3LO_3(_D[j][s-1][s-1][s-1][0], k, _P1nsm, _P2nsm, _P3nsm);
								_D[j][s][s][s][0][k] = fac1 + fac2;
								_D[j][s][s][s][0][k] /= gamma;
							}

							// RecRel #3:
							for (int m=s-1; m>=0; m--)
							{
								double fac1 = (-_alpha_s.Beta1() - r1*_alpha_s.Beta2() - r1*r1*_alpha_s.Beta3()) * _D[j][s][s][m+1][1][k];
								double fac2 = RecRelN3LO_2(_D[j][s-1][s-1][m][0], k, _P1nsm, _P2nsm, _P3nsm);
								_D[j][s][s][m][0][k] = fac1 + fac2;
								_D[j][s][s][m][0][k] /= gamma;
							}

							// RecRel #4:
							for (int t=s-1; t>=0; t--)
							{
								for (int m=t; m>=0; m--) // this one doesn't need to be an int, but removes compiler warning (wrong! (for now))
								{
									double fac0 = -2*b*gamma*_D[j][s][t+1][m+1][1][k];
									double fac1 = fac0 - 2.0*(-b*_alpha_s.Beta1() - r1*_alpha_s.Beta1() + c*_alpha_s.Beta2() + c*r1*_alpha_s.Beta3())*_D[j][s][t+1][m+1][1][k];
									double fac2 = RecRelN3LO_4(_D[j][s-1][t][m][0], k, _P1nsm, _P2nsm, _P3nsm);
									_D[j][s][t][m][0][k] = fac1 + fac2;
									_D[j][s][t][m][0][k] /= gamma;
								}
							}
							
						}

						for (uint j=32; j<31+_nf; j++)
						{
							for (uint t=1; t<=s; t++)
							{
								for (uint m=1; m<=t; m++)
								{
									for (uint n=1; n<=m; n++)
										_D[j][s][t][m][n][k] = RecRelN3LO_1(_D[j][s-1][t-1][m-1][n-1], k, _P0ns);
								}
							}

							{
								double fac1 = 0.5*(_alpha_s.Beta1() + r1*_alpha_s.Beta2() - fac*_alpha_s.Beta3()) *_D[j][s][s][s][1][k];
								double fac2 = RecRelN3LO_3(_D[j][s-1][s-1][s-1][0], k, _P1nsp, _P2nsp, _P3nsp);
								_D[j][s][s][s][0][k] = fac1 + fac2;
								_D[j][s][s][s][0][k] /= gamma;
							}

							for (int m=s-1; m>=0; m--)
							{
								double fac1 = (-_alpha_s.Beta1() - r1*_alpha_s.Beta2() - r1*r1*_alpha_s.Beta3()) * _D[j][s][s][m+1][1][k];
								double fac2 = RecRelN3LO_2(_D[j][s-1][s-1][m][0], k, _P1nsp, _P2nsp, _P3nsp);
								_D[j][s][s][m][0][k] = fac1 + fac2;
								_D[j][s][s][m][0][k] /= gamma;
							}

							for (int t=s-1; t>=0; t--)
							{
								for (int m=t; m>=0; m--)
								{
									double fac0 = -2*b*gamma*_D[j][s][t+1][m+1][1][k];
									double fac1 = fac0 - 2.0*(-b*_alpha_s.Beta1() - r1*_alpha_s.Beta1() + c*_alpha_s.Beta2() + c*r1*_alpha_s.Beta3())*_D[j][s][t+1][m+1][1][k];
									double fac2 = RecRelN3LO_4(_D[j][s-1][t][m][0], k, _P1nsp, _P2nsp, _P3nsp);
									_D[j][s][t][m][0][k] = fac1 + fac2;
									_D[j][s][t][m][0][k] /= gamma;
								}
							}
						}

						{
							for (uint t=1; t<=s; t++)
							{
								for (uint m=1; m<=t; m++)
								{
									for (uint n=1; n<=m; n++)
										_D[25][s][t][m][n][k] = RecRelN3LO_1(_D[25][s-1][t-1][m-1][n-1], k, _P0ns);
								}
							}

							{
								double fac1 = 0.5*(_alpha_s.Beta1() + r1*_alpha_s.Beta2() - fac*_alpha_s.Beta3()) *_D[25][s][s][s][1][k];
								double fac2 = RecRelN3LO_3(_D[25][s-1][s-1][s-1][0], k, _P1nsm, _P2nsv, _P3nsv);
								_D[25][s][s][s][0][k] =  fac1 + fac2;
								_D[25][s][s][s][0][k] /= gamma;
							}

							for (int m=s-1; m>=0; m--)
							{
								double fac1 = (-_alpha_s.Beta1() - r1*_alpha_s.Beta2() - r1*r1*_alpha_s.Beta3()) * _D[25][s][s][m+1][1][k];
								double fac2 = RecRelN3LO_2(_D[25][s-1][s-1][m][0], k, _P1nsm, _P2nsv, _P3nsv);
								_D[25][s][s][m][0][k] =  fac1 + fac2;
								_D[25][s][s][m][0][k] /= gamma;
							}

							for (int t=s-1; t>=0; t--)
							{
								for (int m=t; m>=t; m--)
								{
									double fac0 = -2*b*gamma*_D[25][s][t+1][m+1][1][k];
									double fac1 = fac0 - 2.0*(-b*_alpha_s.Beta1() - r1*_alpha_s.Beta1() + c*_alpha_s.Beta2() + c*r1*_alpha_s.Beta3())*_D[25][s][t+1][m+1][1][k];
									double fac2 = RecRelN3LO_4(_D[25][s-1][t][m][0], k, _P1nsm, _P2nsv, _P3nsv);
									_D[25][s][t][m][0][k] = fac1 + fac2;
									_D[25][s][t][m][0][k] /= gamma;
								}
							}
						}
					}

					tf = Clock::now();
					std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(tf-t0);
					std::cerr << "[DGLAP] EvolveNonSinglet(): N3LO iteration took " << elapsed.count() << " seconds.\n";
				}				
			} break;
			default:
			{
				std::cerr << "[DGLAP] EvolveNonSinglet() Encountered invalid order: " << _order << '\n';
				exit(1);
			}
		}
	}


	void DGLAPSolver::Resum(const double Q)
	{
		std::cerr << "[DGLAP] Resum(): resumming to tabulated energy\n";
		_alpha1 = _alpha_s.Evaluate(_alpha_s.Masses(_nf), Q, _alpha0);

		double L1 = std::log(_alpha1/_alpha0);
		double L2 = 0.0;
		double L3 = 0.0;
		double L4 = 0.0;

		if (_order == 1)
		{
			L2 = std::log((_alpha1*_alpha_s.Beta1() + 4.0*PI*_alpha_s.Beta0())
							  /(_alpha0*_alpha_s.Beta1() + 4.0*PI*_alpha_s.Beta0()));
		}
		else if (_order == 2)
		{
			L2 = std::log((16.0*PI_2*_alpha_s.Beta0()+4.*PI*_alpha1*_alpha_s.Beta1()+_alpha1*_alpha1*_alpha_s.Beta2())
						  /(16.*PI_2*_alpha_s.Beta0()+4.*PI*_alpha0*_alpha_s.Beta1()+_alpha0*_alpha0*_alpha_s.Beta2()));
			
			// analytic continuation for arctan
			double aux=4.0*_alpha_s.Beta0()*_alpha_s.Beta2() - _alpha_s.Beta1()*_alpha_s.Beta1();
			if (aux>=0)
			{
				L3 = std::atan(2.0*PI*(_alpha1-_alpha0)*std::sqrt(aux)/(2.*PI*(8.*PI*_alpha_s.Beta0()+(_alpha1+_alpha0))+_alpha1*_alpha0*_alpha_s.Beta2()))/std::sqrt(aux);
			}
			else
			{
				L3 = std::tanh(2.*PI*(_alpha1-_alpha0)*std::sqrt(-aux)/(2.*PI*(8.*PI*_alpha_s.Beta0()+(_alpha1+_alpha0))+_alpha1*_alpha0*_alpha_s.Beta2()))/std::sqrt(-aux);
			}
		}
		else if (_order == 3)
		{
			// shorthand
			double r1 = _r1[_nf];
			double b = _b[_nf];
			double c = _c[_nf];

			L2 = std::log((_alpha1 - r1)/(_alpha0 - r1));
			L3 = std::log((_alpha1*_alpha1 + b*_alpha1 + c) / (_alpha0*_alpha0 + b*_alpha0 + c));
			double aux = std::sqrt(-b*b + 4*c); // never negative, no need for analytic continuation
			L4 = std::atan((_alpha1-_alpha0)*aux / (2*_alpha0*_alpha1 + (_alpha0+_alpha1)*b + 2.0*c))/aux;
			// L4 = std::atan((_alpha1-_alpha0)*aux / (2*_alpha0*_alpha1 + (_alpha0+_alpha1)*b + 2.0*c));
		}
		
		std::cerr << "[DGLAP] Resum(): alpha0=" << _alpha0 << ", _alpha1 = " << _alpha1 << '\n';
		std::cerr << "[DGLAP] Resum(): L1=" << L1 << ", L2=" << L2 << ", L3=" << L3 << ", L4=" << L4 << '\n';
		// singlet
		for (uint j=0; j<=1; j++)
		{
			for (uint k=0; k<_grid.Size()-1; k++)
			{
				_F[j*31][k] = _S[0][j][0][k];
					
				for (uint n=1; n<=ITERATIONS-1; n++)
				{
					for (uint t=0; t<=TRUNC_IDX; t++)
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
					
					for (uint n=1; n<=ITERATIONS-1; n++)
					{
						for (uint j=13; j<=12+_nf; j++)
							_F[j][k] += _A[j][n][k]*std::pow(L1, n)/Factorial(n);
						for (uint j=32; j<=30+_nf; j++)
							_F[j][k] += _A[j][n][k]*std::pow(L1, n)/Factorial(n);
					}
				}
			} break;
			case 1: // NLO
			{
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					for (uint j=13; j<=30+_nf; j++)
					{
						_F[j][k] = _B[j][0][0][k];
						if (j == (12+_nf))
							j = 31;
					}
					for (uint s=1; s<=ITERATIONS-1; s++)
					{
						for (uint n=0; n<=s; n++)
						{
							for (uint j=13; j<=12+_nf; j++)
							{
								_F[j][k] += _B[j][s][n][k]*std::pow(L1,n)*std::pow(L2,(s-n))
										  /Factorial(n)/Factorial(s-n);
							}

							for (uint j=32; j<=30+_nf; j++)
							{
								_F[j][k] += _B[j][s][n][k]*std::pow(L1,n)*std::pow(L2,(s-n))
										  /Factorial(n)/Factorial(s-n);
							}
						}
					}
					
				}
			} break;
			case 2: // NNLO
			{
				for (uint k=0; k<_grid.Size()-1;k++)
				{
					for (uint j=25; j<=30+_nf; j++)
					{
						_F[j][k] = _C[j][0][0][0][k];
						if (j == (24+_nf))
							j = 31;
					}
					for (uint s=1; s<=ITERATIONS-1; s++)
					{
						for (uint t=0; t<=s; t++)
						{
							for (uint n=0; n<=t; n++)
							{
								for (uint j=25; j<=24+_nf;j++)
								{
									double orig = _C[j][s][t][n][k];
									double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
									double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
									_F[j][k] += orig*powers/factorials;
								}
								for (uint j=32; j<=30+_nf; j++) {
									double orig = _C[j][s][t][n][k];
									double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
									double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
									_F[j][k] += orig*powers/factorials;
								}
							}
						}
					}
				}
			} break;
			case 3:
			{
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					for (uint j=25; j<=30+_nf; j++)
					{
						_F[j][k] = _D[j][0][0][0][0][k];
						if (j == (24+_nf))
							j = 31;
					}
					for (uint s=1; s<=ITERATIONS-1; s++)
					{
						for (uint t=0; t<=s; t++)
						{
							for (uint m=0; m<=t; m++)
							{
								for (uint n=0; n<=m; n++)
								{
									for (uint j=25; j<=24+_nf;j++)
									{
										double orig = _D[j][s][t][m][n][k];
										double powers = std::pow(L1,n)*std::pow(L2,(m-n))*std::pow(L3,(t-m))*std::pow(L4,(s-t));
										double factorials = Factorial(n)*Factorial(m-n)*Factorial(t-m)*Factorial(s-t);
										_F[j][k] += orig*powers/factorials;
									}
									for (uint j=32; j<=30+_nf;j++)
									{
										double orig = _D[j][s][t][m][n][k];
										double powers = std::pow(L1,n)*std::pow(L2,(m-n))*std::pow(L3,(t-m))*std::pow(L4,(s-t));
										double factorials = Factorial(n)*Factorial(m-n)*Factorial(t-m)*Factorial(s-t);
										_F[j][k] += orig*powers/factorials;
									}
								}
							}
						}
					}
				}
			} break;
		}
		// flavor reconstruction
		double Nf = static_cast<double>(_nf);
		for (uint k=0; k<_grid.Size()-1;k++)
		{
			if (_order>=2)
			{
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

			for (uint j=1; j<=_nf; j++)
			{
				_F[j][k]   =0.5*(_F[j+18][k] + _F[j+12][k]);
				_F[j+6][k] =0.5*(_F[j+18][k] - _F[j+12][k]);
			}

			if (_order<2)
			{
				_F[25][k]=0.0;
				for (uint j=13; j<=12+_nf; j++)
					_F[25][k] += _F[j][k];

				for (uint j=26; j<=24+_nf; j++)
					_F[j][k] = _F[13][k] - _F[j-12][k];
			}
		}

		std::ostringstream filepath_ss;
		filepath_ss << "out-";
		switch(_order)
		{
			case 0:  filepath_ss << "LO"; break;
			case 1:	 filepath_ss << "NLO"; break;
			case 2:	 filepath_ss << "NNLO"; break;
			case 3:	 filepath_ss << "N3LO"; break;
			default: filepath_ss << "ERROR"; break;
		}
		filepath_ss << ".dat";
		std::string filepath = filepath_ss.str();

		FILE* outfile = fopen(filepath.c_str(), "w");

		// for (uint k=0; k<_grid.Size(); k++) {
		// 	fprintf(outfile, "%15.9g", _grid.At(k));
		// 	for (uint j=26; j<=30; j++)
		// 		fprintf(outfile, "  %15.9g", _F[j][k]);
		// 	for (uint j=32; j<=36; j++)
		// 		fprintf(outfile, "  %15.9g", _F[j][k]);	
		// 	fprintf(outfile, "\n");
		// }

		for (uint k=0; k<_grid.Size(); k++) {
		fprintf(outfile, "%15.9g", _grid.At(k));
			for (uint j=0; j<=12; j++)
				fprintf(outfile, "  %15.9g", _F[j][k]);
			fprintf(outfile, "\n");
		}

		fclose(outfile);
	}

	void DGLAPSolver::ResumThreshold()
	{
		// now sum the recursive series for threshold Q values i guess?
		_alpha1 = _alpha_s.Pre(_nf+1);
		double Nf = static_cast<double>(_nf);

		double L1 = std::log(_alpha1/_alpha0);
		double L2 = 0.0;
		double L3 = 0.0;
		double L4 = 0.0;

		if (_order == 1)
		{
			L2 = std::log((_alpha1*_alpha_s.Beta1() + 4.0*PI*_alpha_s.Beta0())
						  /(_alpha0*_alpha_s.Beta1() + 4.0*PI*_alpha_s.Beta0()));
		}
		else if (_order == 2)
		{
			L2 = std::log((16.0*PI_2*_alpha_s.Beta0()+4.*PI*_alpha1*_alpha_s.Beta1()+_alpha1*_alpha1*_alpha_s.Beta2())
						  /(16.*PI_2*_alpha_s.Beta0()+4.*PI*_alpha0*_alpha_s.Beta1()+_alpha0*_alpha0*_alpha_s.Beta2()));
				
			// analytic continuation for arctan
			double aux=4.0*_alpha_s.Beta0()*_alpha_s.Beta2() - _alpha_s.Beta1()*_alpha_s.Beta1();
			if (aux>=0)
			{
				L3 = std::atan(2.0*PI*(_alpha1-_alpha0)*std::sqrt(aux)/(2.*PI*(8.*PI*_alpha_s.Beta0()+(_alpha1+_alpha0))+_alpha1*_alpha0*_alpha_s.Beta2()))/std::sqrt(aux);
			}
			else
			{
				L3 = std::tanh(2.*PI*(_alpha1-_alpha0)*std::sqrt(-aux)/(2.*PI*(8.*PI*_alpha_s.Beta0()+(_alpha1+_alpha0))+_alpha1*_alpha0*_alpha_s.Beta2()))/std::sqrt(-aux);
			}
		}
		else if (_order == 3)
		{
			// shorthand
			double r1 = _r1[_nf];
			double b = _b[_nf];
			double c = _c[_nf];

			L2 = std::log((_alpha1 - r1)/(_alpha0 - r1));
			L3 = std::log((_alpha1*_alpha1 + b*_alpha1 + c) / (_alpha0*_alpha0 + b*_alpha0 + c));
			double aux = std::sqrt(-b*b + 4*c); // never negative, no need for analytic continuation
			L4 = std::atan((_alpha1-_alpha0)*aux / (2*_alpha0*_alpha1 + (_alpha0+_alpha1)*b + 2.0*c))/aux;
			// L4 = std::atan((_alpha1-_alpha0)*aux / (2*_alpha0*_alpha1 + (_alpha0+_alpha1)*b + 2.0*c));
		}

		std::cerr << "[DGLAP] ResumThreshold(): alpha0=" << _alpha0 << ", _alpha1 = " << _alpha1 << '\n';
		std::cerr << "[DGLAP] ResumThreshold(): L1=" << L1 << ", L2=" << L2 << ", L3=" << L3 << ", L4=" << L4 << '\n';

		// singlet 
		for (uint j=0; j<=1; j++)
		{
			for (uint k=0; k<_grid.Size()-1; k++)
			{
				for (uint n=1; n<=ITERATIONS-1; n++)
				{
					for (uint t=0; t<=TRUNC_IDX; t++)
						_S[0][j][0][k] += _S[t][j][n][k]*std::pow(_alpha1, t)*std::pow(L1, n)/Factorial(n);
				}
			}
		}


		// non-singlet
		switch (_order)
		{
			case 0: // LO
			{
				// resum
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					for (uint n=1; n<=ITERATIONS-1; n++)
					{
						for (uint j=13; j<=12+_nf; j++)
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
					for (uint s=1; s<=ITERATIONS-1; s++)
					{
						for (uint n=0; n<=s; n++)
						{
							for (uint j=13; j<=12+_nf; j++)
							{
								_B[j][0][0][k] += _B[j][s][n][k]*std::pow(L1,n)*std::pow(L2,(s-n))
												/Factorial(n)/Factorial(s-n);
							}

							for (uint j=32; j<=30+_nf; j++)
							{
								_B[j][0][0][k] += _B[j][s][n][k]*std::pow(L1,n)*std::pow(L2,(s-n))
												/Factorial(n)/Factorial(s-n);
							}
						}
					}
				}

				// set distributions
				for (uint k=0; k<_grid.Size()-1;k++)
				{
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
			} break;
			case 2: // NNLO
			{
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					for (uint s=1; s<=ITERATIONS-1; s++)
					{
						for (uint t=0; t<=s; t++)
						{
							for (uint n=0; n<=t; n++)
							{
								for (uint j=25; j<=24+_nf; j++)
								{
									double orig = _C[j][s][t][n][k];
									double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
									double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
									_C[j][0][0][0][k] += orig*powers/factorials;
								}
								for (uint j=32; j<=30+_nf; j++)
								{
									double orig = _C[j][s][t][n][k];
									double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
									double factorials = Factorial(n)*Factorial(t-n)*Factorial(s-t);
									_C[j][0][0][0][k] += orig*powers/factorials;
								}
							}
						}
					}
				}

				for (uint k=0; k<_grid.Size()-1;k++)
				{
					_C[13][0][0][0][k] = _C[25][0][0][0][k];
					for (uint j=26; j<=24+_nf; j++)
						_C[13][0][0][0][k] += _C[j][0][0][0][k];
					_C[13][0][0][0][k] /= Nf;

					for (uint j=14;j<=12+_nf;j++)
						_C[j][0][0][0][k] = _C[13][0][0][0][k] - _C[j+12][0][0][0][k];

					_C[19][0][0][0][k] = _S[0][1][0][k];
					for (uint j=32; j<=30+_nf; j++)
						_C[19][0][0][0][k] += _C[j][0][0][0][k];
					_C[19][0][0][0][k] /= Nf;

					for (uint j=20; j<=18+_nf; j++)
						_C[j][0][0][0][k] = _C[19][0][0][0][k] - _C[j+12][0][0][0][k];

					for (uint j=1; j<=_nf; j++)
					{
						_C[j][0][0][0][k]  =0.5*(_C[j+18][0][0][0][k]+_C[j+12][0][0][0][k]);
						_C[j+6][0][0][0][k]=0.5*(_C[j+18][0][0][0][k]-_C[j+12][0][0][0][k]);
					}
				}
			} break;
			case 3: // N3LO
			{
				for (uint k=0; k<_grid.Size()-1; k++)
				{
					for (uint s=1; s<ITERATIONS; s++)
					{
						for (uint t=0; t<=s; t++)
						{
							for (uint m=0; m<=t; m++)
							{
								for (uint n=0; n<=m; n++)
								{
									for (uint j=25; j<=24+_nf; j++)
									{
										double orig = _D[j][s][t][m][n][k];
										double powers = std::pow(L1,n)*std::pow(L2,(m-n))*std::pow(L3,(t-m))*std::pow(L4,(s-t));
										double factorials = Factorial(n)*Factorial(m-n)*Factorial(t-m)*Factorial(s-t);
										_D[j][0][0][0][0][k] += orig*powers/factorials;
									}

									for (uint j=32; j<=30+_nf; j++)
									{
										double orig = _D[j][s][t][m][n][k];
										double powers = std::pow(L1,n)*std::pow(L2,(m-n))*std::pow(L3,(t-m))*std::pow(L4,(s-t));
										double factorials = Factorial(n)/Factorial(m-n)/Factorial(t-m)/Factorial(s-t);
										_D[j][0][0][0][0][k] += orig*powers/factorials;
									}
								}
							}
						}
					}
				}
				
				for (uint k=0; k<_grid.Size()-1;k++)
				{
					_D[13][0][0][0][0][k] = _D[25][0][0][0][0][k];
					for (uint j=26; j<=24+_nf; j++)
						_D[13][0][0][0][0][k] += _D[j][0][0][0][0][k];
					_D[13][0][0][0][0][k] /= Nf;

					for (uint j=14;j<=12+_nf;j++)
						_D[j][0][0][0][0][k] = _D[13][0][0][0][0][k] - _D[j+12][0][0][0][0][k];

					_D[19][0][0][0][0][k] = _S[0][1][0][k];
					for (uint j=32; j<=30+_nf; j++)
						_D[19][0][0][0][0][k] += _D[j][0][0][0][0][k];
					_D[19][0][0][0][0][k] /= Nf;

					for (uint j=20; j<=18+_nf; j++)
						_D[j][0][0][0][0][k] = _D[19][0][0][0][0][k] - _D[j+12][0][0][0][0][k];

					for (uint j=1; j<=_nf; j++)
					{
						_D[j][0][0][0][0][k]  =0.5*(_D[j+18][0][0][0][0][k]+_D[j+12][0][0][0][0][k]);
						_D[j+6][0][0][0][0][k]=0.5*(_D[j+18][0][0][0][0][k]-_D[j+12][0][0][0][0][k]);
					}
				}
			}
		}
	}


	void DGLAPSolver::HeavyFlavorTreatment()
	{
		return;
		std::cerr << "[DGLAP] HeavyFlavorTreatment(): " << _nf+1 << "th quark mass threshold (mass "
				  << _alpha_s.Masses(_nf+1) << ")\n";

		//Copy of pre-threshold distributions
		for (uint k=0; k<_grid.Size()-1; k++)
		{
			for (uint j=0; j<=1; j++)
				_S[0][j][1][k] = _S[0][j][0][k];


			for (uint i=1; i<=_nf; i++)
			{
				for (uint j=i; j<=i+6; j+=6)
				{
					if (_order == 2)
						_C[j][1][0][0][k]    = _C[j][0][0][0][k];
					else if (_order == 3)
						_D[j][1][0][0][0][k] = _D[j][0][0][0][0][k];
				}
			}
		}

		double aux = _alpha_s.Post(_nf+1);
		std::cerr << "[DGLAP] HeavyFlavorTreatment(): value of alpha_s post threshold: " << aux << '\n';

		//Computation of after-threshold distributions
		for (uint k=0; k<_grid.Size()-1;k++)
		{
			for (uint i=1; i<=_nf; i++)
			{
				for (uint j=i; j<=i+6; j+=6)
				{
					if (_order == 2)
						_C[j][0][0][0][k]    += std::pow(aux/(4.0*PI), 2) * _grid.Convolution(_C[j][1][0][0], _A2ns, k);
					else if (_order == 3)
						_D[j][0][0][0][0][k] += std::pow(aux/(4.0*PI), 2) * _grid.Convolution(_D[j][1][0][0][0], _A2ns, k);
				}
			}

			_S[0][0][0][k] += std::pow(aux/(4.0*PI), 2) * (_grid.Convolution(_S[0][1][1], _A2gq, k)
				            + _grid.Convolution(_S[0][0][1], _A2gg, k));

			double fac = 0.5*std::pow(aux/(4.0*PI), 2) * (_grid.Convolution(_S[0][1][1], _A2hq, k)
				       + _grid.Convolution(_S[0][0][1], _A2hg, k));

			if (_order == 2)
				_C[_nf+1][0][0][0][k]    = _C[_nf+7][0][0][0][k] = fac;
			else if (_order == 3)
				_D[_nf+1][0][0][0][0][k] = _D[_nf+7][0][0][0][0][k] = fac;
		}
	}
	
	
} // namespace Candia2
