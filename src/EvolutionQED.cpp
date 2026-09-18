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
		auto const& alphaqed = _alpha_qed.value();

		// qcd kernels
		auto& p0ns = getExpression(ExprName::P0ns);
		auto& p0qq = getExpression(ExprName::P0qq);
		auto& p0qg = getExpression(ExprName::P0qg);
		auto& p0gq = getExpression(ExprName::P0gq);
		auto& p0gg = getExpression(ExprName::P0gg);

		// qed kernels
		auto& p0ff = getExpression(ExprName::P0ff);
		auto& p0uu = getExpression(ExprName::P0uu);
		auto& p0dd = getExpression(ExprName::P0dd);
		auto& p0ll = getExpression(ExprName::P0ll);
		auto& p0fy = getExpression(ExprName::P0fy);
		auto& p0yf = getExpression(ExprName::P0yf);
		auto& p0yy = getExpression(ExprName::P0yy);

		auto deltaud = static_cast<uint>(QEDPartonIndices::DELTAUD);
		auto sigma = static_cast<uint>(QEDPartonIndices::SIGMA);
		auto gluon = static_cast<uint>(QEDPartonIndices::G);
		auto photon = static_cast<uint>(QEDPartonIndices::PHOTON);
		auto sigmal = static_cast<uint>(QEDPartonIndices::SIGMAL);
		std::array s_dists{deltaud, sigma, gluon, photon, sigmal};

		// most of these are zero with our initial conditions
		// but I wanted to go ahead and add all of them so that
		// once the initial conditions allow non-vanishing valence lepton distributions
		// i won't forget
		std::array ns_dists{
			std::vector{
				static_cast<uint>(QEDPartonIndices::UV),
				// static_cast<uint>(QEDPartonIndices::CV),
				// static_cast<uint>(QEDPartonIndices::DELTAUC),
			},
			std::vector{
				static_cast<uint>(QEDPartonIndices::DV),
				// static_cast<uint>(QEDPartonIndices::SV),
				// static_cast<uint>(QEDPartonIndices::BV),
				// static_cast<uint>(QEDPartonIndices::DELTADS),
				// static_cast<uint>(QEDPartonIndices::DELTASB),
			},
			std::vector{
				// static_cast<uint>(QEDPartonIndices::EV),
				// static_cast<uint>(QEDPartonIndices::MUV),
				// static_cast<uint>(QEDPartonIndices::TAUV),
				static_cast<uint>(QEDPartonIndices::DELTAL2),
				static_cast<uint>(QEDPartonIndices::DELTAL3),	
			}
		};
		std::array ns_splitfuncs{
			p0uu,p0dd,p0ll
		};
		auto num_ns_splitfunc_options = ns_splitfuncs.size();
		
		for (uint j : s_dists)
			std::ranges::copy(_S_QED(0,j,0), arr_singlet.get()[j].begin());
		for (uint j : ns_dists | std::views::join)
			std::ranges::copy(_A_QED(0,j,0), arr_ns.get()[j].begin());

		// for nf=4, we have an equal number of up/down quarks
		// so the numerator, Nup-Ndown = 0, and deltaNf = 0
		double deltaNf = static_cast<double>(alphaqed.numUp()-alphaqed.numDown())/_nf;
		double beta0qcd = _alpha_s.beta0();
		double beta0qed = alphaqed.beta0();
		constexpr double cp = 0.5*((2./3.)*(2./3.) + (-1./3.)*(-1./3.));
		constexpr double cm = 0.5*((2./3.)*(2./3.) - (-1./3.)*(-1./3.));
		double fac_qcd = -2.0/beta0qcd;
		double fac_qed = -2.0/beta0qed;
		fac_qed = 0.0;
		L0QED = 0;

		log(LOG_DEBUG, "evolveQED()", "{:>10} = {: }", "deltaNf", deltaNf);
		log(LOG_DEBUG, "evolveQED()", "{:>10} = {: }", "beta0QED", beta0qed);
		log(LOG_DEBUG, "evolveQED()", "{:>10} = {: }", "beta0QCD", beta0qcd);
		log(LOG_DEBUG, "evolveQED()", "{:>10} = {: }", "cp", cp);
		log(LOG_DEBUG, "evolveQED()", "{:>10} = {: }", "cm", cm);
		log(LOG_DEBUG, "evolveQED()", "{:>10} = {: }", "fac_qcd", fac_qcd);
		log(LOG_DEBUG, "evolveQED()", "{:>10} = {: }", "fac_qed", fac_qed);

		// singlet
		{
		    for (uint s=1; s<_iterations; s++) {
				for (uint n=1; n<=s; n++) {
					auto pows = std::pow(L0QED,n)*std::pow(L0QCD,s-n)/factorial(n)/factorial(s-n);
					for (uint k=0; k<_grid.size()-1;k++) {
						double res1 = fac_qed*(
							cp*_grid.convolution(_S_QED(deltaud,0,n-1), p0ff, k) +
							cm*_grid.convolution(_S_QED(sigma,0,n-1), p0ff, k) +
							0 +
							2*NC*_nf*(cp*deltaNf + cm)*_grid.convolution(_S_QED(photon,0,n-1), p0fy, k) +
							0
						);
						double res2 = fac_qed*(
							cm*_grid.convolution(_S_QED(deltaud,0,n-1), p0ff, k) +
							cp*_grid.convolution(_S_QED(sigma,0,n-1), p0ff, k) +
							0 +
							2*NC*_nf*(cp + cm*deltaNf)*_grid.convolution(_S_QED(photon,0,n-1), p0fy, k) +
							0
						);
						double res4 = fac_qed*(
							cm*_grid.convolution(_S_QED(deltaud,0,n-1), p0yf, k) +
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
									
						_S_QED(deltaud,1,n,k) = res1;
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
							_grid.convolution(_S_QED(deltaud,0,n), p0qq, k) +
							0 +
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
								
						_S_QED(deltaud,1,n,k) = res1;
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
		
		// non-singlet
		{
		    for (uint i=0; i<num_ns_splitfunc_options; ++i) {
				auto const& dists = ns_dists[i];
				auto& p0ff = ns_splitfuncs[i];

				for (uint j : dists) {
					for (uint s=1; s<_iterations; s++) {
						for (uint n=1; n<=s; n++) {
							double pows = std::pow(L0QED,n)*std::pow(L0QCD,s-n)/factorial(n)/factorial(s-n);
							for (uint k=0; k<_grid.size()-1;k++) {
								_A_QED(j,1,n,k) = fac_qed*_grid.convolution(_A_QED(j,0,n-1), p0ff, k);
								arr_ns.get()[j][k] += _A_QED(j,1,n,k)*pows;
							}
						}

						{
							uint n = 0;
							double pows = std::pow(L0QED,n)*std::pow(L0QCD,s-n)/factorial(n)/factorial(s-n);
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
