#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Math.hpp"

#include <cmath>
#include <thread>
#include <print>

namespace Candia2
{
	void DGLAPSolver::evolveNonSingletThreaded(
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
		double L1, double L2, double L3, double L4)
	{
		switch (_order)
		{
			case 0: // LO
			{
				std::println("[THREAD] Performing LO non-singlet evolution threaded.");

                for (uint j=13; j<=12+_nf; ++j)
                    arr.get()[j] = _A2[j][0];
                for (uint j=32; j<=30+_nf; ++j)
                    arr.get()[j] = _A2[j][0];
                    
				std::vector<std::thread> threads{};

                std::println("=============== BEGIN THREAD OUTPUT ====================");
				for (uint j=13; j<=12+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_LO, this, arr, j, L1);	
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_LO, this, arr, j, L1);

			    for (std::thread & t : threads)
					t.join();

                std::println("=============== END THREAD OUTPUT ====================");
				std::println("[THREAD] Finished performing threaded LO non-singlet evolution.");
			} break;
			case 1: // NLO
			{
				std::println("[THREAD] Performing NLO non-singlet evolution threaded.");

                for (uint j=13; j<=12+_nf; ++j)
                    arr.get()[j] = _B2[j][0][0];
                for (uint j=32; j<=30+_nf; ++j)
                    arr.get()[j] = _B2[j][0][0];

				std::vector<std::thread> threads{};
                std::array<double, 2> L{L1, L2};

                std::println("=============== BEGIN THREAD OUTPUT ====================");
				for (uint j=13; j<=12+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NLO, this, arr, j, "P1nsm", L);
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NLO, this, arr, j, "P1nsp", L);

				for (std::thread & t : threads)
					t.join();

                std::println("=============== END THREAD OUTPUT ====================");
				std::println("[THREAD] Finished performing threaded NLO non-singlet evolution.");
			} break;
			case 2: // NNLO
			{
				std::println("[THREAD] Performing NNLO non-singlet evolution threaded.");
				std::vector<std::thread> threads{};

                std::array<double, 3> L{L1, L2, L3};
				std::array<std::string, 2> nsm{"P1nsm", "P2nsm"};
				std::array<std::string, 2> nsp{"P1nsp", "P2nsp"};
				std::array<std::string, 2> nsv{"P1nsm", "P2nsv"};

                std::println("=============== BEGIN THREAD OUTPUT ====================");
				for (uint j=26; j<=24+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NNLO, this, arr, j, nsm, L);
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NNLO, this, arr, j, nsp, L);
				threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_NNLO, this, arr, 25, nsv, L);
				
				for (std::thread & t : threads)
					t.join();

                std::println("=============== END THREAD OUTPUT ====================");
				std::println("[THREAD] Finished performing threaded NNLO non-singlet evolution.");
			} break;
			case 3: // N3LO nonsinglet
			{
				std::println("[THREAD] Performing N3LO non-singlet evolution threaded.");

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

				std::println("=============== BEGIN THREAD OUTPUT ====================");
				for (uint j=26; j<=24+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_N3LO, this, arr, j, nsm, L);
				for (uint j=32; j<=30+_nf; j++)
					threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_N3LO, this, arr, j, nsp, L);
				threads.emplace_back(&DGLAPSolver::_mt_EvolveDistribution_NS_N3LO, this, arr, 25, nsv, L);
				
				for (std::thread & t : threads)
					t.join();

				std::println("=============== END THREAD OUTPUT ====================");
				std::println("[THREAD] Finished performing threaded N3LO non-singlet evolution.");
			} break;
		}
	}


    void DGLAPSolver::_mt_EvolveDistribution_NS_LO (
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
        uint j, double L1)
    {
        for (uint n=0; n<_iterations-1; n++)
		{
            std::println("  [j={}] Iteration {}/{}", j, n, _iterations-1);
			for (uint k=0; k<_grid.size()-1; k++)
            {
				_A2[j][1][k] = recrelLO(_A2[j][0], k, getSplitFunc("P0ns"));
                arr.get()[j][k] += _A2[j][1][k]*std::pow(L1, n)/factorial(n);
            }
			_A2[j][0] = _A2[j][1];
		}
    }
    void DGLAPSolver::_mt_EvolveDistribution_NS_NLO(
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
        uint j, std::string const& P1, std::array<double, 2> const&  L)
    {
        double const L1 = L[0];
        double const L2 = L[1];
        for (uint s=1; s<_iterations; s++)
        {
            std::println("  [j={}] Iteration {}/{}", j, s, _iterations-1);
            for (uint k=0; k<_grid.size()-1;k++)
            {
                for (uint n=1; n<=s; n++)
                {
                    _B2[j][1][n][k] = recrelNLO_1(_B2[j][0][n-1], k, getSplitFunc("P0ns"));
                    arr.get()[j][k] += _B2[j][1][n][k]*std::pow(L1,n)*std::pow(L2,s-n)/factorial(n)/factorial(s-n);
                }

                uint n = 0;
                double res = recrelNLO_2(_B2[j][0][0], k, 
                    getSplitFunc("P0ns"), 
                    getSplitFunc(P1));
                _B2[j][1][0][k] = -_B2[j][1][1][k] + res;
                arr.get()[j][k] += _B2[j][1][0][k]
                    *std::pow(L1,n)*std::pow(L2,s-n)
                    /factorial(n)/factorial(s-n);
            }
            for (uint n=0; n<=s; ++n)
                _B2[j][0][n] = _B2[j][1][n];
        }
    }
    void DGLAPSolver::_mt_EvolveDistribution_NS_NNLO(
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
        uint j, std::array<std::string, 2> const& P, std::array<double, 3> const& L)
    {
        double const L1 = L[0];
        double const L2 = L[1];
        double const L3 = L[2];

        for (uint s=1; s<_iterations; s++)
        {
            std::println("  [j={}] Iteration {}/{}", j, s, _iterations-1);
            for (uint k=0; k<_grid.size()-1; k++)
            {
                // recrel #1:
                for (uint t=1; t<=s; t++)
                {
                    for (uint n=1; n<=t; n++)
                    {
                        double recrel = recrelNNLO_1(_C2[j][0][t-1][n-1], k, getSplitFunc("P0ns"));
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
                        getSplitFunc("P0ns"), getSplitFunc(P[0]), getSplitFunc(P[1]));
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
                        getSplitFunc("P0ns"), getSplitFunc(P[0]));
                    _C2[j][1][t][0][k] = fac1 + fac2;

                    uint n = 0;
                    double orig = _C2[j][1][t][0][k];
                    double powers = std::pow(L1,n)*std::pow(L2,(t-n))*std::pow(L3,(s-t));
                    double factorials = factorial(n)*factorial(t-n)*factorial(s-t);
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
    void DGLAPSolver::_mt_EvolveDistribution_NS_N3LO(
        std::reference_wrapper<std::vector<ArrayGrid>> arr, 
			uint j, std::array<std::string, 3> const& P, std::array<double, 4> const& L)
    {
		double const L1 = L[0];
        double const L2 = L[1];
        double const L3 = L[2];
        double const L4 = L[3];

        // some shorthand
        double r1 = _r1[_nf];
        double b = _b[_nf];
        double c = _c[_nf];
        double gamma = (r1*r1 + r1*b + c)*_alpha_s.beta3();

        for (uint s=1; s<_iterations; s++)
        {
			std::println("  [j={}] Iteration {}/{}", j, s, _iterations-1);
            for (uint k=0; k<_grid.size()-1; k++)
            {
                // recrel #1:
                for (uint t=1; t<=s; t++)
                {
                    for (uint m=1; m<=t; m++)
                    {
                        for (uint n=1; n<=m; n++)
                        {
                            _D2[j][1][t][m][n][k] = recrelN3LO_1(_D2[j][0][t-1][m-1][n-1], k, getSplitFunc("P0ns"));

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
                            factorial(n)
                            *factorial(m-n)
                            *factorial(t-m)
                            *factorial(s-t);
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
