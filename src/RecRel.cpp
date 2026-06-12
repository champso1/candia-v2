// RecRel.cpp

#include "Candia-v2/Candia.hpp"

namespace Candia2
{

	double DGLAPSolver::shift_p1(Expression& P1, Expression& P0, ArrayGridView A, uint k)
	{
		const double L = _log_muf2_mur2;
		const double b0 = _alpha_s.beta0();
		
		double conv1 = _grid.convolution(A, P1, k);
		double conv2 = _grid.convolution(A, P0, k);
		return
			conv1
			- 0.5*b0*L*conv2;
	}
	double DGLAPSolver::shift_p2(Expression& P2, Expression& P1, Expression& P0, ArrayGridView A, uint k)
	{
		double conv1 = _grid.convolution(A, P2, k);
		double conv2 = _grid.convolution(A, P1, k);
		double conv3 = _grid.convolution(A, P0, k);

		const double L = _log_muf2_mur2;
		const double b0 = _alpha_s.beta0();
		const double b1 = _alpha_s.beta1();
		
		return
			conv1
			- b0*L*conv2
			+ 0.25*(b0*b0*L*L - b1*L)*conv3;
	}
	double DGLAPSolver::shift_p3(Expression& P3, Expression& P2, Expression& P1, Expression& P0, ArrayGridView A, uint k)
	{
		double conv1 = _grid.convolution(A, P3, k);
		double conv2 = _grid.convolution(A, P2, k);
		double conv3 = _grid.convolution(A, P1, k);
		double conv4 = _grid.convolution(A, P0, k);

		const double L = _log_muf2_mur2;
		const double b0 = _alpha_s.beta0();
		const double b1 = _alpha_s.beta1();
		const double b2 = _alpha_s.beta2();
		
		return
			conv1
			- 1.5*b0*L*conv2
			+ (0.75*b0*b0*L*L - 0.5*b1*L)*conv3
			+ 0.125*(-b0*b0*b0*L*L*L + 2.5*b0*b1*L*L - b2*L)*conv4;
	}
	
	double DGLAPSolver::recrelS_1(
		ArrayGridView S,
		uint k,
		Expression& P)
	{
		double conv = _grid.convolution(S, P, k);
		return -(2.0/_alpha_s.beta0()) * conv;
	}
	
	double DGLAPSolver::recrelS_2(
		ArrayGridView S_i,
		ArrayGridView S_im1,
		uint k,
		Expression& P0,
		Expression& P1)
	{
		double conv1 = _grid.convolution(S_i, P0, k);

		double fac1 = 2.0/_alpha_s.beta0();
		double fac2 = fac1/TWOPI;
		
		double res = fac1*conv1;

		if (_is_scale_difference) {
			res += fac2*shift_p1(P1, P0, S_im1, k);
		} else {
			double conv2 = _grid.convolution(S_im1, P1, k);
			res += fac2*conv2;
		}
	    
		return -res;
	}

	double DGLAPSolver::recrelS_3(
		ArrayGridView S_i,   // C
		ArrayGridView S_im1, // B
		ArrayGridView S_im2, // A
		uint k,
		Expression& P0,
		Expression& P1,
		Expression& P2)
	{
		// S1 -> S_i
		// S2 -> S_im1
		// S3 -> S_im2
		double conv1 = _grid.convolution(S_i, P0, k);

		double fac1 = 2.0/_alpha_s.beta0();
		double fac2 = 1.0/(PI*_alpha_s.beta0());
		double fac3 = 1.0/(2.0*PI_2*_alpha_s.beta0());
		
		double res = fac1*conv1;

		if (_is_scale_difference) {
			res += fac2*shift_p1(P1, P0, S_im1, k);
			res += fac3*shift_p2(P2, P1, P0, S_im2, k);
		} else {
			double conv2 = _grid.convolution(S_im1, P1, k);
			double conv3 = _grid.convolution(S_im2, P2, k);
			res += fac2*conv2;
			res += fac3*conv3;
		}		
		
		return -res;
	}

	double DGLAPSolver::recrelS_4(
		ArrayGridView S_i,
		ArrayGridView S_im1,
		ArrayGridView S_im2,
		ArrayGridView S_im3,
		uint k,
		Expression& P0,
		Expression& P1,
		Expression& P2,
		Expression& P3)
	{
		double conv1 = _grid.convolution(S_i, P0, k);

		double fac1 = 2.0/_alpha_s.beta0();
		double fac2 = 1.0/(PI*_alpha_s.beta0());
		double fac3 = 1.0/(2.0*PI_2*_alpha_s.beta0());
		double fac4 = 1.0/(4.0*PI_3*_alpha_s.beta0());

		double res = fac1*conv1;

		if (_is_scale_difference) {
			res += fac2*shift_p1(P1, P0, S_im1, k);
			res += fac3*shift_p2(P2, P1, P0, S_im2, k);
			res += fac4*shift_p3(P3, P2, P1, P0, S_im3, k);
		} else {
			double conv2 = _grid.convolution(S_im1, P1, k);
			double conv3 = _grid.convolution(S_im2, P2, k);
			double conv4 = _grid.convolution(S_im3, P3, k);
			res += fac2*conv2;
			res += fac3*conv3;
			res += fac4*conv4;
		}
		
		return -res;
	}



	double DGLAPSolver::recrelLO(
		ArrayGridView A,
		uint k,
		Expression & P0)
	{
		double conv = _grid.convolution(A, P0, k);
		return (-2.0/_alpha_s.beta0())*conv;
	}


	double DGLAPSolver::recrelNLO_1(
		ArrayGridView B,
		uint k,
		Expression& P0)
	{
		double conv = _grid.convolution(B, P0, k);
		return -(2.0/_alpha_s.beta0())*conv;
	}
	double DGLAPSolver::recrelNLO_2(
		ArrayGridView B,
		uint k,
		Expression& P0,
		Expression& P1)
	{

		double fac = -4.0/_alpha_s.beta1();
		double res = 0.0;

		if (_is_scale_difference) {
			res += fac*shift_p1(P1, P0, B, k);
		} else {
		    double conv = _grid.convolution(B, P1, k);
			res += fac*conv;
		}

		return res;
	}


	double DGLAPSolver::recrelNNLO_1(
		ArrayGridView C,
		uint k,
		Expression& P0)
	{
		double conv = _grid.convolution(C, P0, k);
		return -(2.0/_alpha_s.beta0())*conv;
	}
	double DGLAPSolver::recrelNNLO_2(
		ArrayGridView C,
		uint k,
		Expression& P0,
		Expression& P1,
		Expression& P2)
	{
		double fac = -4.0/_alpha_s.beta2();
		double res = 0.0;

		if (_is_scale_difference) {
			res += fac*shift_p2(P2, P1, P0, C, k);
		} else {
			double conv = _grid.convolution(C, P2, k);
			res += fac*conv;
		}
			
	    return res;
	}
	double DGLAPSolver::recrelNNLO_3(
		ArrayGridView C,
		uint k,
		Expression& P0,
		Expression& P1)
	{
		double fac = -8.0;
		double res = 0.0;

		if (_is_scale_difference) {
		    res += fac*shift_p1(P1, P0, C, k);
		} else {
			double conv = _grid.convolution(C, P1, k);
			res += fac*conv;
		}
		
		return res;
	}



	double DGLAPSolver::recrelN3LO_1(
		ArrayGridView D,
		uint k,
		Expression& P0)
	{
		double conv = _grid.convolution(D, P0, k);
		return -(2.0/_alpha_s.beta0()) * conv;
	}
	double DGLAPSolver::recrelN3LO_2(
		ArrayGridView D,
		uint k,
		Expression& P0,
		Expression& P1,
		Expression& P2,
		Expression& P3)
	{
		const double r1 = _r1[_nf];
		const double b = _b[_nf];
		const double c = _c[_nf];

		double fac1 = 32.0*PI_2;
		double fac2 = 16.0*PI*r1;
		double fac3 = -8*(c + b*r1);

		double res = 0.0;

		if (_is_scale_difference) {
			res += fac1*shift_p1(P1, P0, D, k);
			res += fac2*shift_p2(P2, P1, P0, D, k);
			res += fac3*shift_p3(P3, P2, P1, P0, D, k);
		} else {
			double conv1 = _grid.convolution(D, P1, k);
			double conv2 = _grid.convolution(D, P2, k);
			double conv3 = _grid.convolution(D, P3, k);
			res += fac1*conv1;
			res += fac2*conv2;
			res += fac3*conv3;
		}
		
		return res;
	}
	double DGLAPSolver::recrelN3LO_3(
		ArrayGridView D,
		uint k,
		Expression& P0,
		Expression& P1,
		Expression& P2,
		Expression& P3)
	{
		const double r1 = _r1[_nf];
		[[maybe_unused]] const double b = _b[_nf];
		[[maybe_unused]] const double c = _c[_nf];

		double conv0 =  _grid.convolution(D, P0, k);

		double fac1 = -64*PI_2;
		double fac2 = -32*PI*r1;
		double fac3 = -16*r1*r1;
		
		double res = 0.0;
		
		if (_is_scale_difference) {
			res += fac1*shift_p1(P1, P0, D, k);
			res += fac2*shift_p2(P2, P1, P0, D, k);
			res += fac3*shift_p3(P3, P2, P1, P0, D, k);
		} else {
			double conv1 = _grid.convolution(D, P1, k);
			double conv2 = _grid.convolution(D, P2, k);
			double conv3 = _grid.convolution(D, P3, k);
			res += fac1*conv1;
			res += fac2*conv2;
			res += fac3*conv3;
		}
			
	    return res;
	}
	double DGLAPSolver::recrelN3LO_4(
		ArrayGridView D,
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

		double fac1 = 128*PI_2*(b+r1);
		double fac2 = -64*PI*c;
		double fac3 = -32*c*r1;
		double res = 0.0;

		if (_is_scale_difference) {
			res += fac1*shift_p1(P1, P0, D, k);
			res += fac2*shift_p2(P2, P1, P0, D, k);
			res += fac3*shift_p3(P3, P2, P1, P0, D, k);
		} else {
			double conv1 = _grid.convolution(D, P1, k);
			double conv2 = _grid.convolution(D, P2, k);
			double conv3 = _grid.convolution(D, P3, k);
			res += fac1*conv1;
			res += fac2*conv2;
			res += fac3*conv3;
		}
		
	    return res;
	}

} // namespace Candia2
