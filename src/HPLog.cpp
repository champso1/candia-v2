#include "Candia-v2/HPLog.hpp"
#include "Candia-v2/Common.hpp"

#include <numbers>
#include <iostream>
#include <fstream>
#include <limits>

#include <complex>
using namespace std::complex_literals;

namespace Candia2
{

	HPLog::HPLog(const uint W, std::string const& file_prefix)
		: _W(W), _file_prefix(file_prefix)
	{
		LoadCoeffs();
		LoadWeights();
	}

	void HPLog::LoadCoeffs()
	{
		std::cerr << "[HPLOG] LoadCoeffs(); Reading coefficients....\n";
		_DA2 = ReadFile<double>("TTT2.m", 360);
		_DA3 = ReadFile<double>("TTT3.m", 1600);
		_DA4 = ReadFile<double>("TTT4.m", 5040);
		_DA5 = ReadFile<double>("TTT5.m", 17280);
		_DA6 = ReadFile<double>("TTT6.m", 51040);
		_DA7 = ReadFile<double>("TTT7.m", 162240);
		_DA8 = ReadFile<double>("TTT8.m", 486000);
		std::cerr << "[HPLOG] LoadCoeffs(): Read all coefficients\n";
	}

	void HPLog::LoadWeights()
	{
		std::cerr << "[HPLOG] LoadWeights(): Loading weights...\n";
	
		if (_W < 5)
			return;

		try
		{
			if (_W >= 5)
			{
				_ii15 = ReadFile<int>("B51.m", 48);
				_ii25 = ReadFile<int>("B52.m", 48);
				_ii35 = ReadFile<int>("B53.m", 48);
				_ii45 = ReadFile<int>("B54.m", 48);
				_ii55 = ReadFile<int>("B55.m", 48);
			}
			if (_W >= 6)
			{
				_ii16 = ReadFile<int>("B61.m", 116);
				_ii26 = ReadFile<int>("B62.m", 116);
				_ii36 = ReadFile<int>("B63.m", 116);
				_ii46 = ReadFile<int>("B64.m", 116);
				_ii56 = ReadFile<int>("B65.m", 116);
				_ii66 = ReadFile<int>("B66.m", 116);
			}
			if (_W >= 7)
			{
				_ii17 = ReadFile<int>("B71.m", 312);
				_ii27 = ReadFile<int>("B72.m", 312);
				_ii37 = ReadFile<int>("B73.m", 312);
				_ii47 = ReadFile<int>("B74.m", 312);
				_ii57 = ReadFile<int>("B75.m", 312);
				_ii67 = ReadFile<int>("B76.m", 312);
				_ii77 = ReadFile<int>("B77.m", 312);
			}
			if (_W == 8)
			{
				_ii18 = ReadFile<int>("B81.m", 810);
				_ii28 = ReadFile<int>("B82.m", 810);
				_ii38 = ReadFile<int>("B83.m", 810);
				_ii48 = ReadFile<int>("B84.m", 810);
				_ii58 = ReadFile<int>("B85.m", 810);
				_ii68 = ReadFile<int>("B86.m", 810);
				_ii78 = ReadFile<int>("B87.m", 810);
				_ii88 = ReadFile<int>("B88.m", 810);
			}
		} catch (std::exception& e) {
			std::cerr << e.what() << std::endl;
			exit(1);
		}

		std::cerr << "[HPLOG] LoadWeights(): Weights loaded.\n";
	}


	static void __throw_b_coeff(std::string const& B)
	{
		std::cerr << "[HPLOG] Check_B():" << B << " has invalid coefficients.\n";
		exit(1);
	}

	int HPLog::Check_B2(int i1, int i2) {
		std::vector<int> ii1({0,-1,-1});
		std::vector<int> ii2({1,1,0});

		// index cannot be -1, this is a good error-checking value
		int j = -1;
		for (uint i=0; i<ii1.size(); i++) {
			if (i1 == ii1[i] && i2 == ii2[i])
				j = i;
		}

		if (j == -1)
			__throw_b_coeff("B2");

		return j;
	}


	int HPLog::Check_B3(int i1, int i2, int i3) {
		std::vector<int> ii1({0 ,0 ,-1,-1,-1,-1,-1,-1});
		std::vector<int> ii2({1 ,0 ,1 ,1 ,0 ,0 ,-1,-1});
		std::vector<int> ii3({1 ,1 ,1 ,0 ,1 ,0 ,1 ,0 });

		int j = -1;
		for (uint i=0; i<8; i++) {
			if (i1 == ii1[i] && i2 == ii2[i] && i3 == ii3[i])
				j = i;
		}

		if (j == -1)
			__throw_b_coeff("B3");

		return j;
	}

	int HPLog::Check_B4(int i1, int i2, int i3, int i4) {
		std::vector<int> ii1({0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1});
		std::vector<int> ii2({1,0,0, 1, 1, 1, 1, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1});
		std::vector<int> ii3({1,1,0, 1, 1, 0, 0, 1, 1, 0, 0,-1, 1, 1, 0, 0,-1,-1});
		std::vector<int> ii4({1,1,1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0});

		int j = -1;
		for (uint i=0; i<18; i++) {
			if (i1 == ii1[i] && i2 == ii2[i] && i3 == ii3[i] && i4 == ii4[i])
				j = i;
		}

		if (j == -1)
			__throw_b_coeff("B4");

		return j;
	}

	int HPLog::Check_B5(int i1, int i2, int i3, int i4, int i5) {
		int j = -1;
		for (uint i=0; i<48; i++) {
			if (i1 == _ii15[i] && i2 == _ii25[i] && i3 == _ii35[i] && i4 == _ii45[i] && i5 == _ii55[i])
				j = i;
		}

		if (j == -1)
			__throw_b_coeff("B5");

		return j;
	}

	int HPLog::Check_B6(int i1, int i2, int i3, int i4, int i5, int i6) {
		int j = -1;
		for (uint i=0; i<116; i++) {
			if (i1 == _ii16[i] && i2 == _ii26[i] && i3 == _ii36[i] && i4 == _ii46[i] && i5 == _ii56[i] && i6 == _ii66[i])
				j = i;
		}

		if (j == -1)
			__throw_b_coeff("B6");

		return j;
	}

	int HPLog::Check_B7(int i1, int i2, int i3, int i4, int i5, int i6, int i7) {
		int j = -1;
		for (uint i=0; i<312; i++) {
			if (i1 == _ii17[i] && i2 == _ii27[i] && i3 == _ii37[i] && i4 == _ii47[i] && i5 == _ii57[i] && i6 == _ii67[i] && i7 == _ii77[i])
				j = i;
		}

		if (j == -1)
			__throw_b_coeff("B7");

		return j;
	}

	int HPLog::Check_B8(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8) {
		int j = -1;
		for (uint i=0; i<810; i++) {
			if (i1 == _ii18[i] && i2 == _ii28[i] && i3 == _ii38[i] && i4 == _ii48[i] && i5 == _ii58[i] && i6 == _ii68[i] && i7 == _ii78[i] && i8 == _ii88[i])
				j = i;
		}

		if (j == -1)
			__throw_b_coeff("B8");

		return j;
	}



	uint HPLog::Check_X(double x)
	{
		uint type = 0;
		
		double new_x = x;
		if (new_x < 0.0)
		{
			type |= MappingType::Negative;
			new_x = -new_x;
		}

		if (new_x > 1.0)
		{
			type |= MappingType::GRT1;
			new_x = 1.0/new_x;
		}

		if ((new_x > std::sqrt(2)-1) && (new_x <= 1.0))
		{
			type |= MappingType::GRTSqrt2Minus1;
			new_x = (1.0-x)/(1.0+x);
		}

		// sanity check
		if ((new_x < 0.0) || (new_x > std::sqrt(2.0)+1))
		{
			std::cerr << "[HPLOG] Check_X(): Error determining mappings.\n";
			exit(1);
		}

		return type;
	}



	double HPLog::H1Mappings(uint mappings, int i1, double x)
	{
		switch (mappings)
		{
			case HPLog::MappingType::Negative:
			{
				double y = -x;
				switch (i1)
				{
					case -1:
					{
						return -H1(1,y);
					} break;

					case 0:
					{
						std::complex<double> res = std::numbers::pi*1i + H1(0,y);
						return res.real();
					} break;

					case 1:
					{
						return -H1(-1,y);
					} break;

					default:
					{
						__throw_b_coeff("H1");

						return std::numeric_limits<double>::max();
					}
				}
			} break;
			case HPLog::MappingType::GRT1:
			{
				std::cerr << "[HPLOG] H1(): impossible x-value: " << x << " (abs(x) > 1 not possible)\n";
				exit(1);
			} break;
			case HPLog::MappingType::GRTSqrt2Minus1:
			{
				double y = (1.0-x)/(1.0+x);
				switch (i1)
				{
					case -1:
					{
						return std::numbers::ln2 - H1(-1, y);
					} break;

					case 0:
					{
						return -H1(-1,y) - H1(1,y);
					} break;

					case 1:
					{
						return -std::numbers::ln2 + H1(-1,y) - H1(0,y);
					} break;

					default:
					{
						__throw_b_coeff("H1");

						return std::numeric_limits<double>::max();
					}
				}
			} break;
		}

		// unreachable, but avoids compiler warning
		return std::numeric_limits<double>::max();
	}


	double HPLog::H2Mappings(uint mappings, int i1, int i2, double x)
	{
		switch (mappings)
		{
			case Negative:
			{
				if ((i1 == 0) && (i2 == 1))
				{
				}
				else if ((i1 == -1) && (i2 == 1))
				{
				}
				else if ((i1 == -1) && (i2 == 0))
				{
				}
				else
				{
					std::cerr << "[HPLOG] H2(): i1=" << i1 << " and i2=" << i2 << " are invalid.\n";
					exit(1);
				}
				double y = -x;
				std::complex<double> res = (-std::numbers::pi*1i - H1(0, y))*H1(1, y) + H2(0,1, y);
				return res.real();
			} break;
			case GRT1:
			{
				std::cerr << "[HPLOG] H2(): Impossible x value: " << x << '\n';
				exit(1);
			} break;
			case GRTSqrt2Minus1:
			{
				double y = (1.0-x)/(1.0+x);
				double h = H1(-1, y);
				std::complex<double> res = -Zeta2/2.0 + 0.5*h*h + H2(-1,1, y);
				return res.real();
			} break;
		}

		return std::numeric_limits<double>::max(); 
	}



	double HPLog::H1(int i1, double x) {
		uint mappings = Check_X(x);

		if (mappings != MappingType::None)
			return H1Mappings(mappings, i1, x);

		
		double T = 0.0;

		switch(i1) {
			case -1: T = log(1.0+x); break;
			case 0:  T = log(x);     break;
			case 1:  T = log(1.0-x); break;
			default: __throw_b_coeff("B1/H1");
		}

		return T;
	}


	double HPLog::H2(int i1, int i2, double x) {
		uint mappings = Check_X(x);

		// the only weight-2 polylogarithm we compute is H[-1,0],
		// so we can skip checking for i1 and i2 (but with an assert just in case)
		if ((i1 != -1) || (i2 != 0))
		{
			std::cerr << "[HPLOG] H2(): i1=" << i1 << " and i2=" << i2 << "is not valid.\n";
			exit(1);
		}

		if (mappings != MappingType::None)
			return H2Mappings(mappings, i1, i2, x);
		
		int j = Check_B2(i1, i2);
			
		double u =  std::log1p( x);
		double v = -std::log1p(-x);
		double t1, t2, t3;

		t1 = 0.0;
		for (uint i=0; i<40; i++)
			t1 = t1 + std::pow(u,i+1)*_DA2[j*40 + i];

		t2 = 0.0;
		for (uint i=0; i<40; i++) 
			t2 = t2 + std::pow(u,i+1)*_DA2[120 + j*40 + i];
		t2 *= log(u);
	
		t3 = 0.0;
		for (uint i=0; i<40; i++)
			t3 = t3 + std::pow(v,i+1)*_DA2[240 + j*40+i];

		return t1+t2+t3;
	}



	double HPLog::H3(int i1, int i2, int i3, double x) {
		Check_X(x);

		int j = Check_B3(i1, i2, i3);

		double u =  std::log1p( x);
		double v = -std::log1p(-x);
		double du = log(u);
		double dv = log(v);

		double t1, t2, t3, t4, t5;

		t1 = 0.0;
		for (uint i=0; i<40; i++)
			t1 += std::pow(u,i+1)*_DA3[j*40 + i];

		t2 = 0.0;
		for (uint i=0; i<40; i++)
			t2 += std::pow(u,i+1)*_DA3[320 + j*40 + i];
		t2 *= du;

		t3 = 0.0;
		for (uint i=0; i<40; i++)
			t3 += std::pow(u,i+1)*_DA3[640 + j*40 + i];
		t3 *= du*du;

		t4 = 0.0;
		for (uint i=0; i<40; i++)
			t4 += std::pow(v,i+1)*_DA3[960 + j*40 + i];

		t5 = 0.0;
		for (uint i=0; i<40; i++)
			t5 = t5 + std::pow(v,i+1)*_DA3[1280 + j*40 + i];
		t5 *= dv;

		return t1+t2+t3+t4+t5;
	}


	double HPLog::H4(int i1, int i2, int i3, int i4, double x) {
		Check_X(x);

		int j = Check_B4(i1, i2, i3, i4);

		double u =  std::log1p( x);
		double v = -std::log1p(-x);
		double du = log(u);
		double dv = log(v);

		double t1, t2, t3, t4, t5, t6, t7;

	
		t1=0.0;
		for (uint i=0; i<40; i++)
			t1=t1+std::pow(u,i+1)*_DA4[j*40+i];

		t2=0.0;
		for (uint i=0; i<40; i++)
			t2=t2+std::pow(u,i+1)*_DA4[720+j*40+i];
		t2=t2*du;
	
		t3=0.0;
		for (uint i=0; i<40; i++)
			t3=t3+std::pow(u,i+1)*_DA4[1440+j*40+i];
		t3=t3*du*du;
	
		t4=0.0;
		for (uint i=0; i<40; i++)
			t4=t4+std::pow(u,i+1)*_DA4[2160+j*40+i];
		t4=t4*du*du*du;
	
		t5=0.0;
		for (uint i=0; i<40; i++)
			t5=t5+std::pow(v,i+1)*_DA4[2880+j*40+i];

		t6=0.0;
		for (uint i=0; i<40; i++)
			t6=t6+std::pow(v,i+1)*_DA4[3600+j*40+i];
		t6=t6*dv;
	
		t7=0.0;
		for (uint i=0; i<40; i++)
			t7=t7+std::pow(v,i+1)*_DA4[4320+j*40+i];

		t7=t7*dv*dv;

		return t1+t2+t3+t4+t5+t6+t7;
	}


	double HPLog::H5(int i1, int i2, int i3, int i4, int i5, double x) {
		Check_X(x);

		int j = Check_B5(i1, i2, i3, i4, i5);

		double u =  std::log1p( x);
		double v = -std::log1p(-x);
		double du = log(u);
		double dv = log(v);

		double t1, t2, t3, t4, t5, t6, t7, t8, t9;

		t1=0.0;
		for (uint i=0; i<40; i++)
			t1=t1+std::pow(u,i+1)*_DA5[j*40+i];

		t2=0.0;
		for (uint i=0; i<40; i++)
			t2=t2+std::pow(u,i+1)*_DA5[1920+j*40+i];
		t2=t2*du;
	
		t3=0.0;
		for (uint i=0; i<40; i++)
			t3=t3+std::pow(u,i+1)*_DA5[3840+j*40+i];
		t3=t3*du*du;
	
		t4=0.0;
		for (uint i=0; i<40; i++)
			t4=t4+std::pow(u,i+1)*_DA5[5760+j*40+i];
		t4=t4*du*du*du;
	
		t5=0.0;
		for (uint i=0; i<40; i++)
			t5=t5+std::pow(u,i+1)*_DA5[7680+j*40+i];
		t5=t5*du*du*du*du;
	
		t6=0.0;
		for (uint i=0; i<40; i++)
			t6=t6+std::pow(v,i+1)*_DA5[9600+j*40+i];
	
		t7=0.0;
		for (uint i=0; i<40; i++)
			t7=t7+std::pow(v,i+1)*_DA5[11520+j*40+i];
		t7=t7*dv;

	
		t8=0.0;
		for (uint i=0; i<40; i++)
			t8=t8+std::pow(v,i+1)*_DA5[13440+j*40+i];
		t8=t8*dv*dv;
	
		t9=0.0;
		for (uint i=0; i<40; i++)
			t9=t9+std::pow(v,i+1)*_DA5[15360+j*40+i];
		t9=t9*dv*dv*dv;

		return t1+t2+t3+t4+t5+t6+t7+t8+t9;
	}


	double HPLog::H6(int i1, int i2, int i3, int i4, int i5, int i6, double x) {
		Check_X(x);

		int j = Check_B6(i1, i2, i3, i4, i5, i6);

		double u =  std::log1p( x);
		double v = -std::log1p(-x);
		double du = log(u);
		double dv = log(v);

		double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;

	
		t1 = 0.0;
		for (uint i=0; i<40; i++)
			t1 += std::pow(u,i+1)*_DA6[j*40+i];

		t2=0.0;
		for (uint i=0; i<40; i++)
			t2=t2+std::pow(u,i+1)*_DA6[4640+j*40+i];
		t2=t2*du;
	
		t3=0.0;
		for (uint i=0; i<40; i++)
			t3=t3+std::pow(u,i+1)*_DA6[9280+j*40+i];
		t3=t3*du*du;
	
		t4=0.0;
		for (uint i=0; i<40; i++)
			t4=t4+std::pow(u,i+1)*_DA6[13920+j*40+i];
		t4=t4*du*du*du;
	
		t5=0.0;
		for (uint i=0; i<40; i++)
			t5=t5+std::pow(u,i+1)*_DA6[18560+j*40+i];
		t5=t5*du*du*du*du;
	
		t6=0.0;
		for (uint i=0; i<40; i++)
			t6=t6+std::pow(u,i+1)*_DA6[23200+j*40+i];
		t6=t6*du*du*du*du*du;
	
		t7=0.0;
		for (uint i=0; i<40; i++)
			t7=t7+std::pow(v,i+1)*_DA6[27840+j*40+i];

		t8=0.0;
		for (uint i=0; i<40; i++)
			t8=t8+std::pow(v,i+1)*_DA6[32480+j*40+i];
		t8=t8*dv;
	
		t9=0.0;
		for (uint i=0; i<40; i++)
			t9=t9+std::pow(v,i+1)*_DA6[37120+j*40+i];
		t9=t9*dv*dv;
	
		t10=0.0;
		for (uint i=0; i<40; i++)
			t10=t10+std::pow(v,i+1)*_DA6[41760+j*40+i];
		t10=t10*dv*dv*dv;
	
		t11=0.0;
		for (uint i=0; i<40; i++)
			t11=t11+std::pow(v,i+1)*_DA6[46400+j*40+i];
		t11=t11*dv*dv*dv*dv;

		return t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11;
	}


	double HPLog::H7(int i1, int i2, int i3, int i4, int i5, int i6, int i7, double x) {
		Check_X(x);

		int j = Check_B7(i1, i2, i3, i4, i5, i6, i7);

		double u =  std::log1p( x);
		double v = -std::log1p(-x);
		double du = log(u);
		double dv = log(v);

		double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13;

		t1=0.0;
		for (uint i=0; i<40; i++)
			t1=t1+std::pow(u,i+1)*_DA7[j*40+i];
	
		t2=0.0;
		for (uint i=0; i<40; i++)
			t2=t2+std::pow(u,i+1)*_DA7[12480+j*40+i];
		t2=t2*du;
	
		t3=0.0;
		for (uint i=0; i<40; i++)
			t3=t3+std::pow(u,i+1)*_DA7[24960+j*40+i];
		t3=t3*du*du;
	
		t4=0.0;
		for (uint i=0; i<40; i++)
			t4=t4+std::pow(u,i+1)*_DA7[37440+j*40+i];
		t4=t4*du*du*du;
	
		t5=0.0;
		for (uint i=0; i<40; i++)
			t5=t5+std::pow(u,i+1)*_DA7[49920+j*40+i];
		t5=t5*du*du*du*du;
	
		t6=0.0;
		for (uint i=0; i<40; i++)
			t6=t6+std::pow(u,i+1)*_DA7[62400+j*40+i];
		t6=t6*du*du*du*du*du;
	
		t7=0.0;
		for (uint i=0; i<40; i++)
			t7=t7+std::pow(u,i+1)*_DA7[74880+j*40+i];
		t7=t7*du*du*du*du*du*du;
	
		t8=0.0;
		for (uint i=0; i<40; i++)
			t8=t8+std::pow(v,i+1)*_DA7[87360+j*40+i];

		t9=0.0;
		for (uint i=0; i<40; i++)
			t9=t9+std::pow(v,i+1)*_DA7[99840+j*40+i];
		t9=t9*dv;
	
		t10=0.0;
		for (uint i=0; i<40; i++)
			t10=t10+std::pow(v,i+1)*_DA7[112320+j*40+i];
		t10=t10*dv*dv;
	
		t11=0.0;
		for (uint i=0; i<40; i++)
			t11=t11+std::pow(v,i+1)*_DA7[124800+j*40+i];
		t11=t11*dv*dv*dv;
	
		t12=0.0;
		for (uint i=0; i<40; i++)
			t12=t12+std::pow(v,i+1)*_DA7[137280+j*40+i];
		t12=t12*dv*dv*dv*dv;
	
		t13=0.0;
		for (uint i=0; i<40; i++)
			t13=t13+std::pow(v,i+1)*_DA7[149760+j*40+i];
		t13=t13*dv*dv*dv*dv*dv;

		return t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13;
	}

	double HPLog::H8(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, double x) {
		Check_X(x);

		int j = Check_B8(i1, i2, i3, i4, i5, i6, i7, i8);

		double u =  std::log1p( x);
		double v = -std::log1p(-x);
		double du = log(u);
		double dv = log(v);

		double t1, t2, t3, t4, t5, t6, t7, t8,
			t9, t10, t11, t12, t13, t14, t15;
		std::vector<double> t(16); t[0] = 0.0;

		t1=0.0;
		for (uint i=0; i<40; i++)
			t1=t1+std::pow(u,i+1)*_DA8[j*40+i];
	
		t2=0.0;
		for (uint i=0; i<40; i++)
			t2=t2+std::pow(u,i+1)*_DA8[32400+j*40+i];
		t2=t2*du;
	
		t3=0.0;
		for (uint i=0; i<40; i++)
			t3=t3+std::pow(u,i+1)*_DA8[64800+j*40+i];
		t3=t3*du*du;
	
		t4=0.0;
		for (uint i=0; i<40; i++)
			t4=t4+std::pow(u,i+1)*_DA8[97200+j*40+i];
		t4=t4*du*du*du;
	
		t5=0.0;
		for (uint i=0; i<40; i++)
			t5=t5+std::pow(u,i+1)*_DA8[129600+j*40+i];
		t5=t5*du*du*du*du;
	
		t6=0.0;
		for (uint i=0; i<40; i++)
			t6=t6+std::pow(u,i+1)*_DA8[162000+j*40+i];
		t6=t6*du*du*du*du*du;
	
		t7=0.0;
		for (uint i=0; i<40; i++)
			t7=t7+std::pow(u,i+1)*_DA8[194400+j*40+i];
		t7=t7*du*du*du*du*du*du;

		t8=0.0;
		for (uint i=0; i<40; i++)
			t8=t8+std::pow(u,i+1)*_DA8[226800+j*40+i];
		t8=t8*du*du*du*du*du*du*du;
	
		t9=0.0;
		for (uint i=0; i<40; i++)
			t9=t9+std::pow(v,i+1)*_DA8[259200+j*40+i];
	
		t10=0.0;
		for (uint i=0; i<40; i++)
			t10=t10+std::pow(v,i+1)*_DA8[291600+j*40+i];
		t10=t10*dv;
	
		t11=0.0;
		for (uint i=0; i<40; i++)
			t11=t11+std::pow(v,i+1)*_DA8[324000+j*40+i];
		t11=t11*dv*dv;
	
		t12=0.0;
		for (uint i=0; i<40; i++)
			t12=t12+std::pow(v,i+1)*_DA8[356400+j*40+i];
		t12=t12*dv*dv*dv;
	
		t13=0.0;
		for (uint i=0; i<40; i++)
			t13=t13+std::pow(v,i+1)*_DA8[388800+j*40+i];
		t13=t13*dv*dv*dv*dv;
	
		t14=0.0;
		for (uint i=0; i<40; i++)
			t14=t14+std::pow(v,i+1)*_DA8[421200+j*40+i];
		t14=t14*dv*dv*dv*dv*dv;
	
		t15=0.0;
		for (uint i=0; i<40; i++)
			t15=t15+std::pow(v,i+1)*_DA8[453600+j*40+i];
		t15=t15*dv*dv*dv*dv*dv*dv;

		return t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15;
	}
	
};

