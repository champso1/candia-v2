#ifndef __SPECIALFUNCS_HPP
#define __SPECIALFUNCS_HPP


namespace Candia2
{
	double Li2(double x);
	double Li3(double x);
	double S12(double x);

	// the following functions are implemented in fortran
	// and thus require C ABI linkage
	extern "C"
	{
		struct dcomplex
		{
			double Re, Im;
		};

		double hplog_(
			double *, int *,
			dcomplex *, dcomplex *, dcomplex *, dcomplex *,
			double * ,double *, double *, double *,
			double *, double *, double *, double *,
			int *, int *);
	}
	
}; // namespace Candia2


#endif // __SPECIALFUNCS_HPP
