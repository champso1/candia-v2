# 1.9.0 (from 1.8.0) <date>
## Bugs
	- The Fortran N3LO splitting functions were called with the wrong imod (the default `ImodAvg` selected B), `P3gg` was missing an `else` (so `Imod1` gave the average), and `P3nsv`/`P3gg::calcDelta` ignored imod entirely (was commented out)
	- LHAPDF output and `LHAPDFDistribution` input had u and d (and ubar and dbar) swapped; `LHAPDFDistribution` also copied s into sbar and never filled charm. This was only the special case for the LesHouchesDistribution but was accidentally carried over the other distributions
	-`calculateSubtractionPDFs` subtracted had some mixups between charm/bottom due to some copy-paste mistakes
	- `Grid::interpFindIdx` cached the size of the first grid ever used, so a second grid of different size interpolated incorrectly. This of course didn't affect single runs, so all previous results were valid
   - warnings from `CANDIA_ENABLE_WARNINGS` no longer leak to consumers, `CANDIA_LIBOME_DIR` and `CANDIA_LHAPDF_DIR` work as documented, and `FindLHAPDF` sets `LHAPDF_INCLUDE_DIRS`/`LHAPDF_LIBRARIES` correctly
   - Evolving to the initial energy, or exactly onto a threshold with the upper number of flavours, no longer returns all-zero PDFs.
   - The subtraction PDFs now correctly use the number of light flavours (3 for charm, 4 for bottom) in the OMEs.
   - The subtraction PDFs of a heavy quark are now computed when it is active at the final energy, and left at zero otherwise. The condition was inverted: they were only computed for Q <= m_h, so an evolution to 100 GeV returned none and one to 3 GeV returned the bottom ones but not the charm ones.
   - In the FFNS, whenever muR /= muF, alphas wasn't correctly evaluated at kr*muF=muR.
   
## Non-Major Changes
   - The grid no longer contains x=1 at all, and every loop now covers every grid point. Previous behavior was artifact of candia-v1, where final grid points was determined by final xtab value. This is unecessary, and x=1 is singular anyway, so it's just removed entirely. Another nice effect is that the loops went to size-1 to not bring in the x=1 point, which doesn't happen anymore.
   - Errors via LOG_ERROR now throw a `Candia2::Error` (nothing more than `std::runtime_error`) instead of calling `exit()`, including errors raised in worker threads. They are raised regardless of verbosity.
   - `Candia2::assert` was renamed to `Candia2::ensure`, so that including `<cassert>` before a Candia header compiles.
   - `SplittingFunction::update` and `OpMatElem::update` are now member functions, meaning nf, beta0 and LM are stored per expression rather than shared by every solver in the process, so several solvers can run at the same time.
   - `trunc_idx` below the perturbative order is raised to the order (with a warning) instead of overflowing
   - several input checks and misleading log messages were fixed (e.g. printing "truncated" when doing exact or something like that)
   - `ArrayGrid` copy/move members are defaulted (moves no longer copy).
   - The Gauss-Legendre nodes for the convolutions now span the full mapping range z in [0, 1] instead of [xmin, 1], which skipped the start of every mapping on accident -> minor improvement in benchmark tables
   - `std::execution::par` is now used instead of `par_unseq` since loop bodies are not really vectorizeable in general. no vectorization happened anyway probably, which is corroborated by the unchanged time and outputs, but this is just a bit safer
   
## Major Changes
	- Added CTest regression tests in `tests/`, and all things previously in `tests/` were moved to `studies/` (some of which were deleted, the others still need some major reworks I think). Calling them "tests" (while not having any *actual* tests) always bothered me. this makes more sense, and separates independent/one-off codes from actual tests that should be ran every time changes are made. See README for description on how to build/use.
	- `LHAPDFGrid` now starts every evolution at the PDF's initial energy `Q0()` (previously the first requested Q, which gave an all-nan  grid in the example), splits the grid into subgrids at heavy-quark thresholds inside the Q range (keeping the discontinuity of the PDFs and of alphas), and generates the info file directly from the grid and other input parameters. The infofile template now only needs metadata such as `SetDesc`, `Authors` and `Reference`, and can be omitted.
	- The convolutions with the splitting functions and OMEs are precomputed once per number of flavours as matrices acting on the grid (`ConvolutionKernel`, built by `Grid::buildKernel()`), instead of re-evaluating the expressions and interpolating at every Gauss point for every iteration. Results agree with the previous ones to rounding (about 1e-12 relative to each distribution), and evolutions beyond LO are 20-90x faster: an N3LO evolution to 100 GeV takes about 0.3 s instead of 27 s, and the `lhapdf_grid` example 2 s instead of 95 s.

# 1.8.0 (from 1.7.0) September 24, 2026
	- libome is no longer packaged as a submodule. Previously, this was done because we needed a slight modification to some of its behavior to be compatible with how we do convolutions, but this can get very messy very easily, for instance if the user has their own installation of libome, if there are updates we have to manually patch them in to our fork, etc. So, it is now a separate dependency (i.e. the user must install it themself), and the patch is applied by defining our own "custom" OMEs in LibomeInterface.{hpp,cpp} (whose results have been tested to match the original behavior). We have added `CANDIA_LIBOME_DIR` as a variable to point to the installation prefix.
	- There are non-CANDIA-prefixed versions of the prefix directories for LHAPDF and libome, those being the general for `<lib>_ROOT`, so `libome_ROOT` and `LHAPDF_ROOT` are also available. 
	- Fixed some bugs with missing guards to test if LHAPDF is desired&found before compiling certain files, like the `lhapdf_grid.cpp` example.
	- Cleaned up the FindLHAPDF.cmake file to try and handle the weird case that the recent LHAPDF versions' `lhapdf-config` binaries actually accept the `--includedir` flag instead of `--incdir`, which is what is shown in the help output.
	- All tests now require LHAPDF support to avoid similar issues to what the above point fixed.
	- Removed `yaml-cpp` as a dependency from `tests/main` since it wasn't useful for the end user and wasn't documented or used purposefully.
	- Fixed `lhapdf_grid.cpp` example to use a better LHAPDF function to add the current directory to the PDF set search path.
	- Overall cleanup of the installation interface to, e.g., allow choosing between shared/static library, automatically bringing in LHAPDF if the user added it during compilation, etc.

# 1.7.0 (from 1.6.0) July 11, 2026:
	- Using the exact P3ns splitting functions is now the default, and the `.useP3Exact()` method is removed, with all approximation specifications moving to the `.setP3ApproximationTypes()` method.
	- Specifying to evolve with truncated ansatz is moved from `DGLAPOptions` to just being a new method, `.evolveTrunc()`, with `.evolve()` performing the exact only
	- Fixed issues with the calculation of the subtraction and residual PDFs, as well as the alphas normalization
	- Most LOG\_INFO messages were moved to LOG\_DEBUG, so that by default we only see the evolution progress and nothing else
	- GridFillers were completely removed since LogLinQuad was just objectively better and faster, there was no reason to use anything else. Everything else is moved to `Grid`, and its constructor is now much cleaner
	- Removed `read_table.cpp` from the examples, and replaced it with `lhapdf_grid.cpp`, and cleaned up `evolve_dglap.cpp` to add more comments explaining the different pieces of code

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

