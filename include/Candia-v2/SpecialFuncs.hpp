#ifndef __SPECIALFUNCS_HPP
#define __SPECIALFUNCS_HPP


namespace Candia2
{

	struct dcomplex {
		double Re, Im;
	};
	
	extern double Li2(const double x);
	extern double Li3(const double x);
	extern double S12(const double x);

	// the following functions are implemented in fortran
	extern "C" {
		///! main harmonic polylog calculation routine
		double hplog_(double *, int *,
					  dcomplex *, dcomplex *, dcomplex *, dcomplex *,
					  double * ,double *, double *, double *,
					  double *, double *, double *, double *,
					  int *, int *);

		// NNLO kernels
		double p2nsma_(double*, int*);
		double p2nspa_(double*, int*);
		double p2nssa_(double*, int*);
		double p2nsb_ (double*, int*);
		double p2nsmc_(double*, int*);
		double p2nspc_(double*, int*);
		double p2psa_ (double*, int*);
		double p2qga_ (double*, int*);
		double p2gqa_ (double*, int*);
		double p2gga_ (double*, int*);
		double p2ggb_ (double*, int*);
		double p2ggc_ (double*, int*);
	}
	
}; // namespace Candia2


#endif // __SPECIALFUNCS_HPP
