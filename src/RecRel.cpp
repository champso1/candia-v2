#include "Candia-v2/Candia.hpp"
#include "Candia-v2/FuncArrayGrid.hpp"

namespace Candia2
{	
	double DGLAPSolver::recrelS_1(
		ArrayGrid & S,
		uint k,
		FunctionGrid & P)
	{
		double conv = P.convolution(S, k);
		return -(2.0/_alpha_s.beta0()) * conv;
	}
	
	double DGLAPSolver::recrelS_2(ArrayGrid& S1,
		ArrayGrid& S2,
		uint k,
		FunctionGrid& P0,
		FunctionGrid& P1)
	{
		double conv1 = P0.convolution(S2, k);
		double conv2 = P1.convolution(S1, k);
		
		double res = conv1 * (2.0/_alpha_s.beta0());
		res += conv2 / (PI*_alpha_s.beta0());

		if (_log_mur2_muf2 != 0.0)
		{
		    double convL = P0.convolution(S1, k);
			res -= _log_mur2_muf2*convL/(2.0*PI);
		}
		
		return -res;
	}

	double DGLAPSolver::recrelS_3(ArrayGrid& S1,
		ArrayGrid& S2,
		ArrayGrid& S3,
		uint k,
		FunctionGrid& P0,
		FunctionGrid& P1,
		FunctionGrid& P2)
	{
		double conv1 = P0.convolution(S3, k);
		double conv2 = P1.convolution(S2, k);
		double conv3 = P2.convolution(S1, k);

		double res = conv1 * 2.0/_alpha_s.beta0();
		res += conv2 / (PI*_alpha_s.beta0());
		res += conv3 / (2.0*PI_2*_alpha_s.beta0());

		if (_log_mur2_muf2 != 0.0)
		{
			double convL1 =  P0.convolution(S2, k);
			double convL2 =  P1.convolution(S1, k);
			double convL3 =  P2.convolution(S1, k);

			res -= _log_mur2_muf2*convL1/(2.0*PI);
			res += _log_mur2_muf2*convL2 * (_alpha_s.beta0()*_log_mur2_muf2 - (_alpha_s.beta1()/_alpha_s.beta0()))/(8.0*PI_2);
			res -= _log_mur2_muf2*convL3/(2.0*PI);
		}
		
		
		return -res;
	}

	double DGLAPSolver::recrelS_4(ArrayGrid& S3, // S^{i-3}
		ArrayGrid& S2, // S^{i-2}
		ArrayGrid& S1, // S^{i-1}
		ArrayGrid& S0, // S^{i}
		uint k,
		FunctionGrid& P0,
		FunctionGrid& P1,
		FunctionGrid& P2,
		FunctionGrid& P3)
	{
		double conv1 = P0.convolution(S0, k);
		double conv2 = P1.convolution(S1, k);
		double conv3 = P2.convolution(S2, k);
		double conv4 = P3.convolution(S3, k);

		double res = -conv1 * (2.0/_alpha_s.beta0());
		res -= conv2 / (PI*_alpha_s.beta0());
		res -= conv3 / (2.0*PI_2*_alpha_s.beta0());
		res -= conv4 / (4.0*PI_3*_alpha_s.beta0());

		if (_log_mur2_muf2 != 0.0)
		{
			const double beta0 = _alpha_s.beta0();
			const double beta1 = _alpha_s.beta1();
			const double beta2 = _alpha_s.beta2();
			const double L = _log_mur2_muf2;
			
		    double convL1 =  P0.convolution(S1, k);
			double convL2a = P1.convolution(S2, k);
			double convL2b = P0.convolution(S2, k);
			double convL3a = P2.convolution(S3, k);
			double convL3b = P1.convolution(S3, k);
			double convL3c = P0.convolution(S3, k);

			res += (L/(2.0*PI)) * convL1;
			res += (L/(2.0*PI_2)) * convL2a;
			res += (((beta1/beta0)*L - beta0*L*L)/(8.0*PI_2)) * convL2b;
			
			res += ((3.0*L)/(8.0*PI_3)) * convL3a;
			res += ((2.0*(beta1/beta0)*L - 3.0*beta0*L*L)/(16.0*PI_3)) * convL3b;
			
			res += ((2.0*(beta2/beta0)*L - 5.0*beta1*L*L + 2.0*beta0*beta0*L*L*L)/(64.0*PI_3)) * convL3c;
		}
		
		return res;
	}



	double DGLAPSolver::recrelLO(
		ArrayGrid & A,
		uint k,
		FunctionGrid & P0)
	{
		double conv = P0.convolution(A, k);
		return (-2.0/_alpha_s.beta0())*conv;
	}


	double DGLAPSolver::recrelNLO_1(
		ArrayGrid& B,
		uint k,
		FunctionGrid& P0)
	{
		double conv = P0.convolution(B, k);
		return -(2.0/_alpha_s.beta0())*conv;
	}
	double DGLAPSolver::recrelNLO_2(
		ArrayGrid& B,
		uint k,
		FunctionGrid& P0,
		FunctionGrid& P1)
	{
		double conv = P1.convolution(B, k);
		double res = -(4.0/_alpha_s.beta1())*conv;

		if (_log_mur2_muf2 != 0.0)
		{
			double convL = P0.convolution(B, k);
			res += (2.0*_log_mur2_muf2*_alpha_s.beta0()/_alpha_s.beta1()) * convL;
		}
		
		return res;
	}


	double DGLAPSolver::recrelNNLO_1(
		ArrayGrid& C,
		uint k,
		FunctionGrid& P0)
	{
		double conv = P0.convolution(C, k);
		return -(2.0/_alpha_s.beta0())*conv;
	}
	double DGLAPSolver::recrelNNLO_2(
		ArrayGrid& C,
		uint k,
		FunctionGrid& P0,
		FunctionGrid& P1,
		FunctionGrid& P2)
	{
		double conv = P2.convolution(C, k);
		double res = -4.0/_alpha_s.beta2() * conv;

		if (_log_mur2_muf2 != 0.0)
		{
			const double beta0 = _alpha_s.beta0();
			const double beta1 = _alpha_s.beta1();
			const double beta2 = _alpha_s.beta2();
			
			double convL1 = P1.convolution(C, k);
			double convL2 = P0.convolution(C, k);


			res += (4.0*(beta0/beta2)*_log_mur2_muf2) * convL1;
			res += ((beta1*_log_mur2_muf2 - beta0*beta0*_log_mur2_muf2*_log_mur2_muf2)/beta2) * convL2;
		}
			
	    return res;
	}
	double DGLAPSolver::recrelNNLO_3(
		ArrayGrid& C,
		uint k,
		FunctionGrid& P0,
		FunctionGrid& P1)
	{
		double conv = P1.convolution(C, k);
		double res = -8.0 * conv;

		if (_log_mur2_muf2 != 0.0)
		{
			double convL = P0.convolution(C, k);
			res += (4.0*_alpha_s.beta0()*_log_mur2_muf2) * convL;
		}
		
		return res;
	}



	double DGLAPSolver::recrelN3LO_1(
		ArrayGrid& D,
		uint k,
		FunctionGrid& P0)
	{
		double conv = P0.convolution(D, k);
		return -(2.0/_alpha_s.beta0()) * conv;
	}
	double DGLAPSolver::recrelN3LO_2(
		ArrayGrid& D,
		uint k,
		FunctionGrid& P0,
		FunctionGrid& P1,
		FunctionGrid& P2,
		FunctionGrid& P3)
	{
		const double r1 = _r1[_nf];
		const double b = _b[_nf];
		const double c = _c[_nf];

		double conv0 = P0.convolution(D, k);
		double conv1 = P1.convolution(D, k);
		double conv2 = P2.convolution(D, k);
		double conv3 = P3.convolution(D, k);

		const double fac1 = 32.0*PI_2;
		const double fac2 = 16.0*PI*r1;
		const double fac3 = -8*(c + b*r1);

		double res = fac1*conv1 + fac2*conv2 + fac3*conv3;

		if (_log_mur2_muf2 != 0.0)
		{
			const double beta0 = _alpha_s.beta0();
			const double beta1 = _alpha_s.beta1();
			const double beta2 = _alpha_s.beta2();
			const double L = _log_mur2_muf2;

			res += (-(fac1/2.0)*beta0*L) * conv0;

			res += (-fac2*beta0*L) * conv1;
			res += (-(fac2/4.0)*(beta1*L - beta0*beta0*L*L)) * conv0;

			res += ((-3.0*fac3/2.0)*beta0*L) * conv2;
			res += (-(fac3/4.0)*(2.0*beta1*L - 3.0*beta0*beta0*L*L)) * conv1;
			res += (-(fac3/16.0)*(2.0*beta2*L - 5.0*beta1*beta2*L*L + 2.0*beta0*beta0*beta0*L*L*L)) * conv0;
		}
		
		return res;
	}
	double DGLAPSolver::recrelN3LO_3(
		ArrayGrid& D,
		uint k,
		FunctionGrid& P0,
		FunctionGrid& P1,
		FunctionGrid& P2,
		FunctionGrid& P3)
	{
		const double r1 = _r1[_nf];
		// const double b = _b[_nf];
		// const double c = _c[_nf];

		double conv0 =  P0.convolution(D, k);
		double conv1 =  P1.convolution(D, k);
		double conv2 =  P2.convolution(D, k);
		double conv3 =  P3.convolution(D, k);

		const double fac1 = -64*PI_2;
		const double fac2 = -32*PI*r1;
		const double fac3 = -16*r1*r1;
		double res = fac1*conv1 + fac2*conv2 + fac3*conv3;
		
		if (_log_mur2_muf2 != 0.0)
		{
			const double beta0 = _alpha_s.beta0();
			const double beta1 = _alpha_s.beta1();
			const double beta2 = _alpha_s.beta2();
			const double L = _log_mur2_muf2;

			res += (-(fac1/2.0)*beta0*L) * conv0;

			res += (-fac2*beta0*L) * conv1;
			res += (-(fac2/4.0)*(beta1*L - beta0*beta0*L*L)) * conv0;

			res += ((-3.0*fac3/2.0)*beta0*L) * conv2;
			res += (-(fac3/4.0)*(2.0*beta1*L - 3.0*beta0*beta0*L*L)) * conv1;
			res += (-(fac3/16.0)*(2.0*beta2*L - 5.0*beta1*beta2*L*L + 2.0*beta0*beta0*beta0*L*L*L)) * conv0;
		}
			
	    return res;
	}
	double DGLAPSolver::recrelN3LO_4(
		ArrayGrid& D,
		uint k,
		FunctionGrid& P0,
		FunctionGrid& P1,
		FunctionGrid& P2,
		FunctionGrid& P3)
	{
		const double r1 = _r1[_nf];
		const double b = _b[_nf];
		const double c = _c[_nf];

		double conv0 = P0.convolution(D, k);
		double conv1 = P1.convolution(D, k);
		double conv2 = P2.convolution(D, k);
		double conv3 = P3.convolution(D, k);

		const double fac1 = 128*PI_2*(b+r1);
		const double fac2 = -64*PI*c;
		const double fac3 = -32*c*r1;
		double res = fac1*conv1 + fac2*conv2 + fac3*conv3;

		if (_log_mur2_muf2 != 0.0)
		{
			const double beta0 = _alpha_s.beta0();
			const double beta1 = _alpha_s.beta1();
			const double beta2 = _alpha_s.beta2();
			const double L = _log_mur2_muf2;

			res += (-(fac1/2.0)*beta0*L) * conv0;

			res += (-fac2*beta0*L) * conv1;
			res += (-(fac2/4.0)*(beta1*L - beta0*beta0*L*L)) * conv0;

			res += ((-3.0*fac3/2.0)*beta0*L) * conv2;
			res += (-(fac3/4.0)*(2.0*beta1*L - 3.0*beta0*beta0*L*L)) * conv1;
			res += (-(fac3/16.0)*(2.0*beta2*L - 5.0*beta1*beta2*L*L + 2.0*beta0*beta0*beta0*L*L*L)) * conv0;
		}
		
	    return res;
	}

} // namespace Candia2
