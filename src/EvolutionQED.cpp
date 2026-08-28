// EvolutionThread.cpp

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/ArrayGrid.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/Math.hpp"

#include <cmath>
#include <functional>

namespace Candia2
{
	void DGLAPSolver::evolveQED(
		std::reference_wrapper<std::vector<ArrayGrid>> arr_singlet,
		std::reference_wrapper<std::vector<ArrayGrid>> arr_ns,
		double L0QED, double L0QCD)
	{
		log(LOG_DEBUG, "DGLAP", "Performing the evolution with QED effects");

		auto& p0ns = getExpression(ExprName::P0ns);
		auto& p0qq = getExpression(ExprName::P0qq);
		auto& p0qg = getExpression(ExprName::P0qg);
		auto& p0gq = getExpression(ExprName::P0gq);
		auto& p0gg = getExpression(ExprName::P0gg);
		
		auto& p0uu = getExpression(ExprName::P0uu);
		auto& p0dd = getExpression(ExprName::P0dd);
		auto& p0ll = getExpression(ExprName::P0ll);
		
		auto& p0uy = getExpression(ExprName::P0uy);
		auto& p0dy = getExpression(ExprName::P0dy);
		auto& p0ly = getExpression(ExprName::P0ly);
		
		auto& p0yu = getExpression(ExprName::P0yu);
		auto& p0yd = getExpression(ExprName::P0yd);
		auto& p0yl = getExpression(ExprName::P0yl);
		
		auto& p0yy = getExpression(ExprName::P0yy);

		auto beta0qcd = _alpha_s.beta0();
		auto beta0qed = _alpha_qed.value().beta0();

		// singlet
		/*
		{
			auto sigmaud = static_cast<uint>(QEDPartonIndices::SIGMAUD);
			auto sigma = static_cast<uint>(QEDPartonIndices::SIGMA);
			auto gluon = static_cast<uint>(QEDPartonIndices::G);
			auto photon = static_cast<uint>(QEDPartonIndices::PHOTON);
			auto sigmal = static_cast<uint>(QEDPartonIndices::SIGMAL);
			std::array s_dists{sigmaud, sigma, gluon, photon, sigmal};

			// for nf=4, we have an equal number of up/down quarks
			// so the numerator, Nup-Ndown = 0, and deltaNf = 0
			double deltaNf = 0;

			double cp = 0.5*((2./3.)*(2./3.) + (-1./3.)*(-1./3.));
			double cm = 0.5*((2./3.)*(2./3.) - (-1./3.)*(-1./3.));
			double fac_qcd = -2.0/beta0qcd();
			double fac_qed = -2.0/beta0qed;

			for (uint s=1; s<_iterations; s++) {
				for (uint n=1; n<=s; n++) {
					auto pows = std::pow(L0QED,n)*std::pow(L0QCD,s-n)/factorial(n)/factorial(s-n);
					for (uint k=0; k<_grid.size()-1;k++) {
						double res1 = fac_qed*(
							cp*_grid.convolution(_S_QED(sigmaud,0,n-1), p0ff, k) +
							cm*_grid.convolution(_S_QED(sigma,0,n-1), p0ff, k) +
							0 +
							2*NC*_nf*(cp*deltaNf + cm)*_grid.convolution(_S_QED(photon,0,n-1), p0fy, k) +
							0
						);
						double res2 = fac_qed*(
							cm*_grid.convolution(_S_QED(sigmaud,0,n-1), p0ff, k) +
							cp*_grid.convolution(_S_QED(sigma,0,n-1), p0ff, k) +
							0 +
							2*NC*_nf*(cp + cm*deltaNf)*_grid.convolution(_S_QED(photon,0,n-1), p0fy, k) +
							0
						);
						double res4 = fac_qed*(
							cm*_grid.convolution(_S_QED(sigmaud,0,n-1), p0yf, k) +
							cp*_grid.convolution(_S_QED(sigma,0,n-1), p0yf, k) +
							0 +
							-3.0/4.0*beta0qed*_grid.convolution(_S_QED(photon,0,n-1), p0yy, k) +
							_grid.convolution(_S_QED(sigmal,0,n-1), p0yf, k)
						);
						double res5 = fac_qed*(
							0 +
							0 +
							0 +
							2.0*_nl*_grid.convolution(_S_QED(photon,0,n-1), p0fy, k) +
							_grid.convolution(_S_QED(sigmal,0,n-1), p0ff, k)
						);
									
						_S_QED(sigmaud,1,n,k) = res1;
						_S_QED(sigma,1,n,k) = res2;
						_S_QED(gluon,1,n,k) = 0; // obviously
						_S_QED(photon,1,n,k) = res4;
						_S_QED(sigmal,1,n,k) = res5;

						for (uint j : s_dists)
							arr_singlet.get()[j][k] += _S_QED(j,1,n,k)*pows;
					}
				}

				{
					uint n = 0;
					auto pows = std::pow(L0QED,n)*std::pow(L0QCD,s-n)/factorial(n)/factorial(s-n);
					for (uint k=0; k<_grid.size()-1;k++) {
						double res1 = fac_qcd*(
							_grid.convolution(_S_QED(sigmaud,0,n), p0ns, k) +
							deltaNf*(
								_grid.convolution(_S_QED(sigma,0,n), p0qq, k) +
								_grid.convolution(_S_QED(sigma,0,n), p0ns, k)) +
							deltaNf*_grid.convolution(_S_QED(gluon,0,n), p0qg, k) +
							0 +
							0
						);
						double res2 = fac_qcd*(
							0 +
							_grid.convolution(_S_QED(sigma,0,n), p0qq, k) +
							_grid.convolution(_S_QED(gluon,0,n), p0qg, k) +
							0 +
							0
						);
						double res3 = fac_qcd*(
							0 +
							_grid.convolution(_S_QED(sigma,0,n), p0gq, k) +
							_grid.convolution(_S_QED(gluon,0,n), p0gg, k) +
							0 +
							0
						);
								
						_S_QED(sigmaud,1,n,k) = res1;
						_S_QED(sigma,1,n,k) = res2;
						_S_QED(gluon,1,n,k) = res3;
						_S_QED(photon,1,n,k) = 0; // obviously
						_S_QED(sigmal,1,n,k) = 0; // obviously
								
						for (uint j : s_dists)
							arr_singlet.get()[j][k] += _S_QED(j,1,n,k)*pows;
					}
				}

				for (uint j : s_dists) {
					for (uint n=0; n<=s; ++n)
						std::ranges::copy(_S_QED(j,1,n), _S_QED(j,0,n).begin());
				}
			}
		}
		*/
		
		// non-singlet
		{
			std::array ns_dists{
				std::vector{
					static_cast<uint>(QEDPartonIndices::UV),
					static_cast<uint>(QEDPartonIndices::DELTAUC),
				},
				std::vector{
					static_cast<uint>(QEDPartonIndices::DV),
					static_cast<uint>(QEDPartonIndices::DELTADS),	
					static_cast<uint>(QEDPartonIndices::DELTASB)
				},
			};
			std::array splitfuncs{
				p0uu,p0dd
			};

			for (uint i=0; i<2; ++i) {
				auto const& dists = ns_dists[i];
				auto& p0ff = splitfuncs[i];

				for (uint j : dists) {
					double fac_qcd = -2.0/beta0qcd;
					double fac_qed = -2.0/beta0qed;
			
					for (uint s=1; s<_iterations; s++) {
						for (uint n=1; n<=s; n++) {
							double pows = std::pow(L0QCD,n)*std::pow(L0QED,s-n)/factorial(n)/factorial(s-n);
							for (uint k=0; k<_grid.size()-1;k++) {
								_A_QED(j,1,n,k) = fac_qed*_grid.convolution(_A_QED(j,0,n-1), p0ff, k);
								arr_ns.get()[j][k] += _A_QED(j,1,n,k)*pows;
							}
						}

						{
							uint n = 0;
							double pows = std::pow(L0QCD,n)*std::pow(L0QED,s-n)/factorial(n)/factorial(s-n);
							for (uint k=0; k<_grid.size()-1;k++) {
								_A_QED(j,1,n,k) = fac_qcd*_grid.convolution(_A_QED(j,0,n), p0ns, k);
								arr_ns.get()[j][k] += _A_QED(j,1,n,k)*pows;
							}
						}
				
						for (uint n=0; n<=s; ++n)
							std::ranges::copy(_A_QED(j,1,n), _A_QED(j,0,n).begin());
					}
				}
			}
		}
	}
}
