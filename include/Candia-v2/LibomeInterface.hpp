/**
 *  @file LibomeInferface.hpp
 *  @brief Contains the interface with libome, particularly a modification of how the plus-distribution coefficient is handled
 */

#include "ome/ome_type_aliases.h"
#include "ome/rpd_distribution.h"

#include "ome/AggQ.h"
#include "ome/AgqQ.h"
#include "ome/AQg.h"
#include "ome/AqgQ.h"
#include "ome/AQqPS.h"
#include "ome/AQqPSs.h"
#include "ome/AqqQNSEven.h"
#include "ome/AqqQNSOdd.h"
#include "ome/AqqQPS.h"

namespace ome
{
	namespace Candia2
	{
		template<typename Tnum, typename Tfunc, typename... Trest>
		class func_candia_plusfunc_omx
		{
		public:
			/// Type alias for the numerical type template parameter
			using numeric_type = Tnum;
			/// Type alias for the callable type template parameter
			using element_type = Tfunc;
			/// Boolean type alias indicating that this class has an eval_plus_int method
			using has_eval_plus_int = std::true_type;

			func_candia_plusfunc_omx()
				: target_function_(element_type()) {};

			explicit
			func_candia_plusfunc_omx(element_type target_function)
				: target_function_(target_function) {};

			// this is the updated piece, no logs and no 1/(1-x),
			// since we evaluate this at x=1
		    numeric_type operator()([[maybe_unused]] numeric_type x, Trest... rest) const
			{
				return(target_function_(0.0, rest...));
			};

			numeric_type eval_plus_int(numeric_type x, Trest... rest) const
			{
				using std::log;
				using std::pow;
				int min_power = target_function_.min_power();
				assert((min_power >= 0, "Only non-negative exponents are supported"));

				numeric_type log_omx = log(static_cast<numeric_type>(1)-x);

				// Calculate the integral over the plus function
				size_t num_coeffs = (target_function_.max_power()+1) - min_power;
				std::vector<numeric_type> plusfunc_ints(num_coeffs,
					static_cast<numeric_type>(0));
				// Compute the exponents that appear in the integral over the plus
				// functions
				std::iota(
					plusfunc_ints.begin(),
					plusfunc_ints.end(),
					static_cast<numeric_type>(min_power+1)
				);
				// Evaluate the integrals over the plus functions
				std::transform(
					plusfunc_ints.begin(),
					plusfunc_ints.end(),
					plusfunc_ints.begin(),
					[log_omx](numeric_type e) { return(-pow(log_omx,e)/e); }
				);

				// Combine the integrals over the individual plus functions with the
				// coefficients
				return(target_function_.eval_subst(plusfunc_ints));
			};

		private:
			element_type target_function_;
		};

		template<typename Tnum>
		using candia_plusfunc_coeff = laurent_polynomial<Tnum, Tnum>;
		template<typename Tnum>
		using candia_plusfunc = func_candia_plusfunc_omx<Tnum, candia_plusfunc_coeff<Tnum>>;
		template<typename Tnum>
		using candia_nf_plus = laurent_polynomial<Tnum, candia_plusfunc<Tnum>, Tnum>;
		template<typename Tnum>
		using candia_logm_plus = laurent_polynomial<Tnum, candia_nf_plus<Tnum>, Tnum, Tnum>;
		template<typename Tnum>
		using candia_as_plus = laurent_polynomial<Tnum, candia_logm_plus<Tnum>, Tnum, Tnum, Tnum>;
		template<typename Tnum>
		using candia_as_plus_view = laurent_polynomial_view<Tnum, candia_logm_plus<Tnum>, Tnum, Tnum, Tnum>;

		using candia_ome_type =
			rpd_distribution<
				ome_as_view<double>,
				candia_as_plus_view<double>,
				ome_as_const_view<double>>;

		extern const candia_as_plus<double> candia_AggQ_plus;
		extern const candia_as_plus<double> candia_AqqQNSEven_plus;
		extern const candia_as_plus<double> candia_AqqQNSOdd_plus;
		
	    extern const candia_ome_type AggQ;
		extern const candia_ome_type AgqQ;
		extern const candia_ome_type AQg;
		extern const candia_ome_type AqgQ;
		extern const candia_ome_type AQqPS;
	    extern const candia_ome_type AQqPSs;
		extern const candia_ome_type AqqQNSEven;
	    extern const candia_ome_type AqqQNSOdd;
	    extern const candia_ome_type AqqQPS;
	}
}


