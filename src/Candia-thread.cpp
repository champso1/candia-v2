#include "Candia-v2/Candia.hpp"

#include <cmath>
#include <format>
#include <iostream>
#include <memory>
#include <thread>

namespace Candia2
{
	void DGLAPSolver::EvolveNonSingletThreaded()
	{
		switch (_order)
		{
			case 0: // LO
			{
				std::cout << "[THREAD] !Performing LO non-singlet evolution threaded!\n";
				std::vector<std::thread> threads{};
				
				for (uint j=13; j<=12+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_LO, this, j);	
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_LO, this, j);

			    for (std::thread & t : threads)
					t.join();

				std::cout << "[THREAD] !Finished performing threaded LO non-singlet evolution\n";
			} break;
			case 1: // NLO
			{
				std::cout << "[THREAD] !Performing NLO non-singlet evolution threaded!\n";
				std::vector<std::thread> threads{};

				for (uint j=13; j<=12+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NLO, this, j, _P1nsm);
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NLO, this, j, _P1nsp);

				for (std::thread & t : threads)
					t.join();

				std::cout << "[THREAD] !Finished performing threaded NLO non-singlet evolution\n";
			} break;
			case 2: // NNLO
			{
				std::cout << "[THREAD] !Performing NNLO non-singlet evolution threaded!\n";
				std::vector<std::thread> threads{};

				std::array<std::shared_ptr<Expression>, 2> nsm{_P1nsm, _P2nsm};
				std::array<std::shared_ptr<Expression>, 2> nsp{_P1nsp, _P2nsp};
				std::array<std::shared_ptr<Expression>, 2> nsv{_P1nsm, _P2nsv};

				for (uint j=26; j<=25+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NNLO, this, j, nsm);
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NNLO, this, j, nsp);
				threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NNLO, this, 25, nsv);
				
				for (std::thread & t : threads)
					t.join();

				std::cout << "[THREAD] !Finished performing threaded NNLO non-singlet evolution\n";
			} break;
			case 3: // N3LO nonsinglet
			{
				std::cout << "[THREAD] !Performing N3LO non-singlet evolution threaded!\n";
				std::vector<std::thread> threads{};

				std::array<std::shared_ptr<Expression>, 3> nsm{_P1nsm, _P2nsm, _P3nsm};
				std::array<std::shared_ptr<Expression>, 3> nsp{_P1nsp, _P2nsp, _P3nsp};
				std::array<std::shared_ptr<Expression>, 3> nsv{_P1nsm, _P2nsv, _P3nsv};

				for (uint j=26; j<=25+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_N3LO, this, j, nsm);
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_N3LO, this, j, nsp);
				threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_N3LO, this, 25, nsv);
				
				for (std::thread & t : threads)
					t.join();

				std::cout << "[THREAD] !Finished performing threaded N3LO non-singlet evolution\n";
			} break;
		}
	}

	static bool found_nan = false;

	void DGLAPSolver::_mt_EvolveDistribution_NS_LO(uint j)
	{
		// solve for all A_n(x) coefficients
		for (uint n=0; n<_iterations-1; n++)
		{
			// std::cerr << "[DGLAP] EvolveNonSinglet(): LO Iteration " << n << '\n';

			for (uint k=0; k<_grid.Size()-1; k++)
			{
				double fac = (-2.0/_alpha_s.Beta0());
				double conv = _grid.Convolution(_A[j][n], _P0ns, k);
				double res = fac*conv;
				
				if (std::isnan(res) && !found_nan)
				{
					// construct the output string in this way because
					// `std::cerr << x << y;` isn't guaranteed to execute in one operation
					// so another thread could spit something out between this thread's x and y
					std::string outstr1 = std::format(
						"[DGLAP] _mt_EvolveDistribution_NS_LO(j={}): NAN at n={} and k={}\n", j, n, k
					);
					std::string outstr2 = std::format(
						"    _A[{}][{}][{}]={}, _grid[{}]={}\n",
						j, n, k, _A[j][n][k], k, _grid[k]
					);
					std::string outstr = outstr1+outstr2;
					std::cerr << outstr;
					found_nan = true;
				}
				_A[j][n+1][k] = res;
			}
			
		}
		std::cout << std::format("[DGLAP] ! Distribution j={} finished! \n", j);
	}

	void DGLAPSolver::_mt_EvolveDistribution_NS_NLO(uint j, std::shared_ptr<Expression> P1)
	{							
	    for (uint s=1; s<=_iterations-1; ++s)
		{
			for (uint k=0; k<_grid.Size()-1; ++k)
			{
			    for (uint n=1; n<=s; n++)
				{
					double conv = _grid.Convolution(_B[j][s-1][n-1], _P0ns, k);
					double fac = -(2.0/_alpha_s.Beta0());
					double res = fac*conv;
					_B[j][s][n][k] = res;
				}

				double res1 = -_B[j][s][1][k];
				double conv = _grid.Convolution(_B[j][s-1][0], P1, k);
				double fac = -(4.0/_alpha_s.Beta1());
				double res2 = fac*conv;
			    
				_B[j][s][0][k] = res1 + res2;
			}
		}
		std::cout << std::format("[DGLAP] ! Distribution j={} finished! \n", j);
	}

	void DGLAPSolver::_mt_EvolveDistribution_NS_NNLO(uint j, std::array<std::shared_ptr<Expression>, 2> const& P)
	{
		for (uint s=1; s<=_iterations-1; s++)
		{
			for (uint k=0; k<_grid.Size()-1; k++)
			{
				for (uint t=1; t<=s; t++)
				{
					for (uint n=1; n<=t; n++)
					{
						double conv = _grid.Convolution(_C[j][s-1][t-1][n-1], _P0ns, k);
						double res = -(2.0/_alpha_s.Beta0())*conv;
									
						_C[j][s][t][n][k] = res;
									
					}
				}

				// RecRel #2:
				{
					double fac1 = -0.5*_C[j][s][s][1][k];
					double conv = _grid.Convolution(_C[j][s-1][s-1][0], P[1], k);
					double fac2 = -4.0/_alpha_s.Beta2() * conv;
								
					_C[j][s][s][0][k] = fac1 + fac2;
				}

				// these must be regular ints;
				// unsigned ints, when they are 0 and get --,
				// underflow back to positive 4b,
				// remaining positive and the loop continues (very bad!)

				// RecRel #3:
				for (int t=s-1; t>=0; t--)
				{
					double fac1 = -2.0*_alpha_s.Beta1()*(_C[j][s][t+1][0][k] + _C[j][s][t+1][1][k]);
					double conv = _grid.Convolution(_C[j][s-1][t][0], P[0], k);
					double fac2 = -8.0*conv;
								
					_C[j][s][t][0][k] = fac1 + fac2;
				}
			}
		}
		std::cout << std::format("[DGLAP] ! Distribution j={} finished! \n", j);
	}


	void DGLAPSolver::_mt_EvolveDistribution_NS_N3LO(uint j, std::array<std::shared_ptr<Expression>, 3> const& P)
	{
		double r1 = _r1[_nf];
		double b = _b[_nf];
		double c = _c[_nf];
		double gamma = (r1*r1 + r1*b + c)*_alpha_s.Beta3();
		
		for (uint s=1; s<=_iterations-1; s++)
		{
			for (uint k=0; k<_grid.Size()-1; k++)
			{
				for (uint t=1; t<=s; t++)
				{
					for (uint m=1; m<=t; m++)
					{
						for (uint n=1; n<=m; n++)
						{
							double conv = _grid.Convolution(_D[j][s-1][t-1][m-1][n-1], _P0ns, k);
							double res = -(2.0/_alpha_s.Beta0()) * conv;
							_D[j][s][t][m][n][k] = res;
						}
					}
				}

				// RecRel #2:
				{
					double fac1 = 0.5*(16.0*PI_2*_alpha_s.Beta1() + 4*PI*r1*_alpha_s.Beta2() - (c + b*r1)*_alpha_s.Beta3());

					double conv1 = _grid.Convolution(_D[j][s-1][s-1][s-1][0], P[0], k);
					double conv2 = _grid.Convolution(_D[j][s-1][s-1][s-1][0], P[1], k);
					double conv3 = _grid.Convolution(_D[j][s-1][s-1][s-1][0], P[2], k);

					const double fac2a = 32.0*PI_2*conv1;
					const double fac2b = 16.0*PI*r1*conv2;
					const double fac2c = -8*(c + b*r1)*conv3;
					double fac2 = fac2a+fac2b+fac2c;

					_D[j][s][s][s][0][k] = fac1 * _D[j][s][s][s][1][k] + fac2;
					_D[j][s][s][s][0][k] /= gamma;

				}

				// RecRel #3:
				for (int m=s-1; m>=0; m--)
				{
					double fac1 = -(16.0*PI_2*_alpha_s.Beta1() + 4.0*PI*r1*_alpha_s.Beta2() + r1*r1*_alpha_s.Beta3());

					double conv1 = _grid.Convolution(_D[j][s-1][s-1][m][0], P[0], k);
					double conv2 = _grid.Convolution(_D[j][s-1][s-1][m][0], P[1], k);
					double conv3 = _grid.Convolution(_D[j][s-1][s-1][m][0], P[2], k);

					const double fac2a = -64*PI_2*conv1;
					const double fac2b = -32*PI*r1*conv2;
					const double fac2c = -16*r1*r1*conv3;
					double fac2 = fac2a+fac2b+fac2c;

					_D[j][s][s][m][0][k] = fac1 * _D[j][s][s][m+1][1][k] + fac2;
					_D[j][s][s][m][0][k] /= gamma;
				}

				// RecRel #4:
				for (int t=s-1; t>=0; t--)
				{
					for (int m=t; m>=0; m--)
					{
						double fac1a = -2.0*b*gamma;
						double fac1b = 32*PI_2*(b+r1)*_alpha_s.Beta1() - 8*PI*c*_alpha_s.Beta2() - 2*c*r1*_alpha_s.Beta3();

						double conv1 = _grid.Convolution(_D[j][s-1][t][m][0], P[0], k);
						double conv2 = _grid.Convolution(_D[j][s-1][t][m][0], P[1], k);
						double conv3 = _grid.Convolution(_D[j][s-1][t][m][0], P[2], k);

						const double fac2a = 128*PI_2*(b+r1)*conv1;
						const double fac2b = -64*PI*c*conv2;
						const double fac2c = -32*c*r1*conv3;
						double fac2 = fac2a+fac2b+fac2c;
									
						_D[j][s][t][m][0][k] = fac1a*_D[j][s][t+1][m+1][0][k]
							+ fac1b*_D[j][s][t+1][m+1][1][k] 
							+ fac2;
						_D[j][s][t][m][0][k] /= gamma;
					}
				}
			}
		}
		std::cout << std::format("[DGLAP] ! Distribution j={} finished! \n", j);
	}

	
}
