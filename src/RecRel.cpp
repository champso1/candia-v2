// RecRel.cpp

#include "Candia-v2/Candia.hpp"

namespace Candia2
{	
	double DGLAPSolver::recrelS_1(
		ArrayGrid& S,
		uint k,
		Expression& P)
	{
		double conv = _grid.convolution(S, P, k);
		return -(2.0/_alpha_s.beta0()) * conv;
	}
	
	double DGLAPSolver::recrelS_2(
		ArrayGrid& S_i,
		ArrayGrid& S_im1,
		uint k,
		Expression& P0,
		Expression& P1)
	{
		double conv1 = _grid.convolution(S_i, P0, k);
		double conv2 = _grid.convolution(S_im1, P1, k);
		
		double res = conv1 * (2.0/_alpha_s.beta0());
		res += conv2 / (PI*_alpha_s.beta0());

		if (_is_scale_difference) {
		    double convL = _grid.convolution(S_im1, P0, k);
			res -= _log_muf2_mur2*convL/(2.0*PI);
		}
		
		return -res;
	}

	double DGLAPSolver::recrelS_3(
		ArrayGrid& S_i, // C
		ArrayGrid& S_im1, // B
		ArrayGrid& S_im2, // A
		uint k,
		Expression& P0,
		Expression& P1,
		Expression& P2)
	{
		// S1 -> S_i
		// S2 -> S_im1
		// S3 -> S_im2
		double conv1 = _grid.convolution(S_i, P0, k);
		double conv2 = _grid.convolution(S_im1, P1, k);
		double conv3 = _grid.convolution(S_im2, P2, k);

		double res = conv1 * (2.0/_alpha_s.beta0());
		res += conv2 / (PI*_alpha_s.beta0());
		res += conv3 / (2.0*PI_2*_alpha_s.beta0());

		if (_is_scale_difference) {
			double L = _log_muf2_mur2;
			double beta1 = _alpha_s.beta1();
			double beta0 = _alpha_s.beta0();

			double convL1 =  _grid.convolution(S_im1, P0, k);
			res -= L*convL1/(2.0*PI);

			double convL2a =  _grid.convolution(S_im2, P1, k);
			double convL2b =  _grid.convolution(S_im2, P0, k);
			res -= L*convL2a/(2.0*PI_2);
			res += convL2b * (beta0*L*L - (beta1/beta0)*L)/(8.0*PI_2);
		}
		
		
		return -res;
	}

	double DGLAPSolver::recrelS_4(
		ArrayGrid& S_i,
		ArrayGrid& S_im1,
		ArrayGrid& S_im2,
		ArrayGrid& S_im3,
		uint k,
		Expression& P0,
		Expression& P1,
		Expression& P2,
		Expression& P3)
	{
		double conv1 = _grid.convolution(S_i, P0, k);
		double conv2 = _grid.convolution(S_im1, P1, k);
		double conv3 = _grid.convolution(S_im2, P2, k);
		double conv4 = _grid.convolution(S_im3, P3, k);

		double beta0 = _alpha_s.beta0();
		double fac1 = 2.0/beta0;
		double fac2 = 1.0/(PI*beta0);
		double fac3 = 1.0/(2.0*PI_2*beta0);
		double fac4 = 1.0/(4.0*PI_3*beta0);

		double res = conv1 * fac1;
		res += conv2 * fac2;
		res += conv3 * fac3;
		res += conv4 * fac4;

		if (_is_scale_difference) {
			const double beta0 = _alpha_s.beta0();
			const double beta1 = _alpha_s.beta1();
			const double beta2 = _alpha_s.beta2();
			const double L = _log_muf2_mur2;
			const double L2 = L*L;
			const double L3 = L2*L;
			
		    double convL1 =  _grid.convolution(S_im1, P0, k);
			res -= L*convL1/(2.0*PI);

			double convL2a = _grid.convolution(S_im2, P1, k);
			double convL2b = _grid.convolution(S_im2, P0, k);
			res -= L*convL2a/(2.0*PI_2);
			res += convL2b * (beta0*L2 - (beta1/beta0)*L)/(8.0*PI_2);

			double convL3a = _grid.convolution(S_im3, P2, k);
			double convL3b = _grid.convolution(S_im3, P1, k);
			double convL3c = _grid.convolution(S_im3, P0, k);
			res -= (3.0*L)/(8.0*PI_3) * convL3a;
			res += (3.0*beta0*L2 - (beta1/beta0)*L)/(16.0*PI_3) * convL3b;
			res += (-beta0*beta0*L3 + (5.0/2.0)*beta1*L2 - (beta2/beta0)*L)/(32.0*PI_3) * convL3c;
		}
		
		return -res;
	}



	double DGLAPSolver::recrelLO(
		ArrayGrid & A,
		uint k,
		Expression & P0)
	{
		double conv = _grid.convolution(A, P0, k);
		return (-2.0/_alpha_s.beta0())*conv;
	}


	double DGLAPSolver::recrelNLO_1(
		ArrayGrid& B,
		uint k,
		Expression& P0)
	{
		double conv = _grid.convolution(B, P0, k);
		return -(2.0/_alpha_s.beta0())*conv;
	}
	double DGLAPSolver::recrelNLO_2(
		ArrayGrid& B,
		uint k,
		Expression& P0,
		Expression& P1)
	{
		double conv = _grid.convolution(B, P1, k);
		double res = -(4.0/_alpha_s.beta1())*conv;

		if (_is_scale_difference) {
			double convL = _grid.convolution(B, P0, k);
			res += (2.0*_log_muf2_mur2*_alpha_s.beta0()/_alpha_s.beta1()) * convL;
		}
		
		return res;
	}


	double DGLAPSolver::recrelNNLO_1(
		ArrayGrid& C,
		uint k,
		Expression& P0)
	{
		double conv = _grid.convolution(C, P0, k);
		return -(2.0/_alpha_s.beta0())*conv;
	}
	double DGLAPSolver::recrelNNLO_2(
		ArrayGrid& C,
		uint k,
		Expression& P0,
		Expression& P1,
		Expression& P2)
	{
		double conv = _grid.convolution(C, P2, k);
		double res = -4.0/_alpha_s.beta2() * conv;

		
		if (_is_scale_difference) {
			const double L = _log_muf2_mur2;
			const double beta0 = _alpha_s.beta0();
			const double beta1 = _alpha_s.beta1();
			const double beta2 = _alpha_s.beta2();
			
			double convL1 = _grid.convolution(C, P1, k);
			double convL2 = _grid.convolution(C, P0, k);

			res += L*4.0*(beta0/beta2) * convL1;
			res += ((beta1*L - beta0*beta0*L*L)/beta2) * convL2;
		}
			
	    return res;
	}
	double DGLAPSolver::recrelNNLO_3(
		ArrayGrid& C,
		uint k,
		Expression& P0,
		Expression& P1)
	{
		double conv = _grid.convolution(C, P1, k);
		double res = -8.0 * conv;

		if (_is_scale_difference) {
			double convL = _grid.convolution(C, P0, k);
			res += _log_muf2_mur2*4.0*_alpha_s.beta0() * convL;
		}
		
		return res;
	}



	double DGLAPSolver::recrelN3LO_1(
		ArrayGrid& D,
		uint k,
		Expression& P0)
	{
		double conv = _grid.convolution(D, P0, k);
		return -(2.0/_alpha_s.beta0()) * conv;
	}
	double DGLAPSolver::recrelN3LO_2(
		ArrayGrid& D,
		uint k,
		Expression& P0,
		Expression& P1,
		Expression& P2,
		Expression& P3)
	{
		const double r1 = _r1[_nf];
		const double b = _b[_nf];
		const double c = _c[_nf];

		double conv0 = _grid.convolution(D, P0, k);
		double conv1 = _grid.convolution(D, P1, k);
		double conv2 = _grid.convolution(D, P2, k);
		double conv3 = _grid.convolution(D, P3, k);

		const double fac1 = 32.0*PI_2;
		const double fac2 = 16.0*PI*r1;
		const double fac3 = -8*(c + b*r1);

		double res = fac1*conv1 + fac2*conv2 + fac3*conv3;

		if (_is_scale_difference) {
			const double beta0 = _alpha_s.beta0();
			const double beta1 = _alpha_s.beta1();
			const double beta2 = _alpha_s.beta2();
			const double L = _log_mur2_muf2;

			res -= (fac1*beta0*L/2.0) * conv0;

			res -= (fac2*beta0*L) * conv1;
			res += ((fac2/4.0)*(beta0*beta0*L*L - beta1*L)) * conv0;

			res -= ((3.0/2.0)*fac3*beta0*L) * conv2;
			res += (-fac3*((3.0/4.0)*beta0*beta0*L*L - (1.0/2.0)*beta1*L)) * conv1;
			res += ((fac3/8.0)*(-beta0*beta0*beta0*L*L*L + (5.0/2.0)*beta1*beta2*L*L - beta2*L)) * conv0;
		}
		
		return res;
	}
	double DGLAPSolver::recrelN3LO_3(
		ArrayGrid& D,
		uint k,
		Expression& P0,
		Expression& P1,
		Expression& P2,
		Expression& P3)
	{
		const double r1 = _r1[_nf];
		// const double b = _b[_nf];
		// const double c = _c[_nf];

		double conv0 =  _grid.convolution(D, P0, k);
		double conv1 =  _grid.convolution(D, P1, k);
		double conv2 =  _grid.convolution(D, P2, k);
		double conv3 =  _grid.convolution(D, P3, k);

		const double fac1 = -64*PI_2;
		const double fac2 = -32*PI*r1;
		const double fac3 = -16*r1*r1;
		double res = fac1*conv1 + fac2*conv2 + fac3*conv3;
		
		if (_is_scale_difference) {
			const double beta0 = _alpha_s.beta0();
			const double beta1 = _alpha_s.beta1();
			const double beta2 = _alpha_s.beta2();
			const double L = _log_mur2_muf2;

			res -= (fac1*beta0*L/2.0) * conv0;

			res -= (fac2*beta0*L) * conv1;
			res += ((fac2/4.0)*(beta0*beta0*L*L - beta1*L)) * conv0;

			res -= ((3.0/2.0)*fac3*beta0*L) * conv2;
			res += (-fac3*((3.0/4.0)*beta0*beta0*L*L - (1.0/2.0)*beta1*L)) * conv1;
			res += ((fac3/8.0)*(-beta0*beta0*beta0*L*L*L + (5.0/2.0)*beta1*beta2*L*L - beta2*L)) * conv0;
		}
			
	    return res;
	}
	double DGLAPSolver::recrelN3LO_4(
		ArrayGrid& D,
		uint k,
		Expression& P0,
		Expression& P1,
		Expression& P2,
		Expression& P3)
	{
		const double r1 = _r1[_nf];
		const double b = _b[_nf];
		const double c = _c[_nf];

		double conv0 = _grid.convolution(D, P0, k);
		double conv1 = _grid.convolution(D, P1, k);
		double conv2 = _grid.convolution(D, P2, k);
		double conv3 = _grid.convolution(D, P3, k);

		const double fac1 = 128*PI_2*(b+r1);
		const double fac2 = -64*PI*c;
		const double fac3 = -32*c*r1;
		double res = fac1*conv1 + fac2*conv2 + fac3*conv3;

		if (_is_scale_difference) {
			const double beta0 = _alpha_s.beta0();
			const double beta1 = _alpha_s.beta1();
			const double beta2 = _alpha_s.beta2();
			const double L = _log_mur2_muf2;

			res -= (fac1*beta0*L/2.0) * conv0;

			res -= (fac2*beta0*L) * conv1;
			res += ((fac2/4.0)*(beta0*beta0*L*L - beta1*L)) * conv0;

			res -= ((3.0/2.0)*fac3*beta0*L) * conv2;
			res += (-fac3*((3.0/4.0)*beta0*beta0*L*L - (1.0/2.0)*beta1*L)) * conv1;
			res += ((fac3/8.0)*(-beta0*beta0*beta0*L*L*L + (5.0/2.0)*beta1*beta2*L*L - beta2*L)) * conv0;
		}
		
	    return res;
	}

} // namespace Candia2
