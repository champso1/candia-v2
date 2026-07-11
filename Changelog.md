# 1.7.0 (from 1.6.0) <date>:
	- Using the exact P3ns splitting functions is now the default, and the `.useP3Exact()` method is removed, with all approximation specifications moving to the `.setP3ApproximationTypes()` method.
	- Specifying to evolve with truncated ansatz is moved from `DGLAPOptions` to just being a new method, `.evolveTrunc()`, with `.evolve()` performing the exact only
	- Fixed issues with the calculation of the subtraction and residual PDFs, as well as the alphas normalization
	- Most LOG\_INFO messages were moved to LOG\_DEBUG, so that by default we only see the evolution progress and nothing else
	- GridFillers were completely removed since LogLinQuad was just objectively better and faster, there was no reason to use anything else. Everything else is moved to `Grid`, and its constructor is now much cleaner

# 1.6.0 (from 1.5.2) May 4, 2026:
	- Threading is now the default and there is no option to change this by specifying the `-DENABLE_THREADING=OFF` flag anymore.
	- The ArrayGrid class is much more sophisticated. It should still behave the same to the end user, but now uses only one backing `std::vector` vs. nesting them. This reduces memory fragmentation and consequently page-faults/cache-misses, making the code faster. Access is done with `operator()` for multi-dimensional ArrayGrids, for 1D ones the `operator[]` works exactly as one would expect. There is also the `view()` method which returns an `std::span` to ease in passing data around.
	- All CMake options are now prefixed with `CANDIA_`.
	- Intel Threading Building Blocks is now a dependency to allow usage of C++ execution policies.
	- Based on the way we perform the convolutions, the x-dependence for both the delta and plus pieces of a generic expression are useless. So, we just precalculate, for a given nf, all of the plus/delta values for all expressions and simply grab them when computing the convolution. All nf-dependent quantities in the regular pieces are also precomputed to reduce time spent recalculating stuff.
	- Added the ability to spit out an LHAPDF grid via the `LHAPDFGrid` class. The funtionality is still very primitive, but it works.
	- The API for the Grid class is a little cleaner, and we have removed the possibility to specify to use GSL routines for the convolution/interpolation. The user is also able to relatively easily specify their own method to fill the grid and perform convolutions by making a class that inherits from `GridFillerBase`.
	- The AlphaS code related to the calculation of the coupling before/after quark mass thresholds is not only cleaner but less arbitrary, making it a little safer for distributions that aren't the hard-coded benchmarking distributions.
	- Logging options are now a little more clear, and a single `verbosity` variable handles what is printed. Anything of a lower priority than what is specified via this variable is not printed (thread output is the same as LOG_INFO)
	- Fixed (hopefully all) bugs for compilation/installation/usage across MacOS and Linux.

