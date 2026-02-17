/**
 *  @file SpecialFuncs.hpp
 *  @brief Contains API to access harmonic polylogarithms available in accompanying Fortran code (using the C ABI)
 */

#ifndef __SPECIALFUNCS_HPP
#define __SPECIALFUNCS_HPP

namespace Candia2
{
	/**
	 *  Returns \f$\mathrm{Li}_2(x) = -\int_0^x \frac{\ln(1-y)}{y}\mathrm{d}y\f$
	 */
	double Li2(double x);
	/**
	 *  Returns \f$\mathrm{Li}_3(x) = \int_0^x \frac{\mathrm{Li}_2(y)}{y}\mathrm{d}y\f$
	 */
	double Li3(double x);
	/**
	 *  Returns \f$\mathrm{S}_{1,2}(x) = \frac{1}{2}\int_0^x \frac{\ln^2(1-y)}{y}\mathrm{d}y\f$
	 */
	double S12(double x);

	// the following functions are implemented in fortran
	// and thus require C ABI linkage
	extern "C"
	{
		struct dcomplex
		{
			double Re, Im;
		};

		/**
		 *  @brief evaluates a harmonic polylogarithm ujp to weight 5
		 *  @param x argument of the 1d harmonic polylog (1dHPL)
		 *  @param nw the maximum requested weight of the 1dHPL
		 *  @param Hc1,Hc2,Hc3,Hc4,Hc5 the complex values of the 1dHPL
		 *  @param Hr1,Hr2,Hr3,Hr4,Hr5 the real parts of @a Hc1,Hc2,Hc3,Hc4,Hr5
		 *  @param Hi1,Hi2,Hi3,Hi4,Hi5 the imag parts of @a Hc1,Hc2,Hc3,Hc4,Hi5
		 *  @param n1,n2 The required range of inidices.
		 *         The allowed ranges are (0,1), (-1,0), (-1,1)
		 */
		double hplog5_(
			double* x, int* nw,
			dcomplex* Hc1, dcomplex* Hc2, dcomplex* Hc3, dcomplex* Hc4, dcomplex* Hc5,
			double* Hr1, double* Hr2, double* Hr3, double* Hr4, double* Hr5,
			double* Hi1, double* Hi2, double* Hi3, double* Hi4, double* Hi5,
			int* n1, int* n2);

		double p3gga(double* y, int* nf, int* imod);
		double p3ggb(double* y, int* nf, int* imod);
		double p3ggc(double* y, int* nf, int* imod);

		double p3gga_2410(double* y, int* nf, int* imod);
		double p3ggb_2410(double* y, int* nf, int* imod);
		double p3ggc_2410(double* y, int* nf, int* imod);

		double p3gqa(double* y, int* nf, int* imod);
		double p3gqa_2404(double* y, int* nf, int* imod);
		double p3gqa_2512(double* y, int* nf, int* imod);

		double p3nsma(double* y, int* nf, int* imod);
		double p3nsmb(double* y, int* nf, int* imod);
		double p3nsmc(double* y, int* nf, int* imod);

		double p3nspa(double* y, int* nf, int* imod);
		double p3nspb(double* y, int* nf, int* imod);
		double p3nspc(double* y, int* nf, int* imod);

		double p3nssa(double* y, int* nf, int* imod);
		double p3psa(double* y, int* nf, int* imod);

		double p3qga(double* y, int* nf, int* imod);		
	}
	
}; // namespace Candia2


#endif // __SPECIALFUNCS_HPP
