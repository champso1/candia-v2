#include "Candia-v2/SpecialFuncs.hpp"


namespace Candia2
{

	double Li2(double x) {
		int nw = 2,
			n1 = -1,
			n2 = 1;

		dcomplex HC1[3], HC2[3][3], HC3[3][3][3], HC4[3][3][3][3];
		double HR1[3], HR2[3][3], HR3[3][3][3], HR4[3][3][3][3];
		double HI1[3], HI2[3][3], HI3[3][3][3], HI4[3][3][3][3];

		hplog_(&x, &nw,
			   HC1, *HC2, **HC3, ***HC4,
			   HR1, *HR2, **HR3, ***HR4,
			   HI1, *HI2, **HI3, ***HI4,
			   &n1, &n2);

		return HR2[2][1];
	}

    
	double Li3(double x) {
		int nw = 3,
			n1 = -1,
			n2 = 1;
		dcomplex HC1[3], HC2[3][3], HC3[3][3][3], HC4[3][3][3][3];
		double   HR1[3], HR2[3][3], HR3[3][3][3], HR4[3][3][3][3];
		double   HI1[3], HI2[3][3], HI3[3][3][3], HI4[3][3][3][3];

		hplog_(&x, &nw,
			   HC1, *HC2, **HC3, ***HC4,
			   HR1, *HR2, **HR3, ***HR4,
			   HI1, *HI2, **HI3, ***HI4,
			   &n1, &n2);

		return HR3[2][1][1];
	}


	double S12(double x)
	{
		int nw = 3,
			n1 = -1,
			n2 = 1;
		dcomplex HC1[3], HC2[3][3], HC3[3][3][3], HC4[3][3][3][3];
		double   HR1[3], HR2[3][3], HR3[3][3][3], HR4[3][3][3][3];
		double   HI1[3], HI2[3][3], HI3[3][3][3], HI4[3][3][3][3];

		hplog_(&x, &nw,
			   HC1, *HC2, **HC3, ***HC4,
			   HR1, *HR2, **HR3, ***HR4,
			   HI1, *HI2, **HI3, ***HI4,
			   &n1, &n2);

		return HR3[2][2][1];
	}



}

