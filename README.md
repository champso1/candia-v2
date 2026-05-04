# Candia-v2


`Candia-v2` is a new and improved version of `Candia` (cfr. C. Hampson, M. Guzzi, arXiv:2512.22667; A. Cafarella, M. Guzzi, C. Coriano', Comp.Phys.Comm. 179 2008; A. Cafarella, M. Guzzi, C. Coriano', Nucl.Phys.B748 2006), a computer code to numerically solve DGLAP evolution for collinear PDFs in the x-space up to next-to-next-to-next-to leading order (N3LO) accuracy in perturbative QCD. `Candia-v2` currently uses an approximate version of the 4-loop splitting functions as the calculation of their exact analytical form is still in progress. `Candia` was originally written in C and was only capable of evolution up to next-to-next-to (NNLO) accuracy in QCD. New information on the splitting functions and operator matrix elements as well as advancements to C++ led to the development of `Candia-v2`, which also brings significant optimizations.

## Building

### Prerequisites

- `libome`: fast interface to OMEs. included as a submodule to this repository
- GSL: GNU Scientific Library
- `tbb`: Intel Threading Building Blocks -- allows usage of C++ execution policies to parallelize/vectorize C++ algorithm loops
- Compiler with C++20 support see [here](#compiler-and-system-support) for info about compiler/system choice
- LHAPDF (optional): for interface to/from LHAPDF grids.
- CMake and a build tool e.g. `Make` or `Ninja`
- Doxygen (optional): for building the code documentation
- LaTeX tools (optional): for building the manuscripts

#### Compiler and System Support

Windows is in general not supported except for Windows Subsystem for Linux. However, Intel provides a C/C++ and Fortran compiler and library suite which, coupled with [this platform-independent CMake-based GSL wrapper](https://github.com/ampl/GSL), allows one to use `Candia-v2` on Windows, just without the LHAPDF support. As far as we are aware, LHAPDF doesn't have a widely available fork/wrapper to port the Autotools-based system to something like CMake. Further, LHAPDF uses some POSIX syscalls which would require some manual intervention.

MacOS is supported but LLVM's CLang or GCC are required, because as far as we are aware, Apple's CLang's `libc++` doesn't implement C++ execution policies. As a note, it was a bit challenging to link with LLVM's `libc++`, so we recommend just using GCC on MacOS. Further, it is best to compile GSL and LHAPDF (if using the LHAPDF interface) using the same GCC version to avoid weird compile errors with `Candia-v2`.

### Compiling

Compiling follows the standard CMake procedure:

```bash
mkdir build
cd build
cmake .. <options>
cmake --build .
# optional: cmake --install .
```

Among standard CMake options like `-DCMAKE_BUILD_TYPE` or `-DCMAKE_INSTALL_PREFIX` are the following `Candia-v2` specific options (all of which are prefixed with `CANDIA_`):
- `CANDIA_WITH_LHAPDF` (bool; default: OFF): compile LHAPDF support. 
- `CANDIA_LHAPDF_DIR` (string; default: '/usr/local'): if LHAPDF is installed in a nonstandard location (i.e. not /usr/local), then specify its installation prefix here. this variable is ignored if `CANDIA_WITH_LHAPDF` is OFF.
- `CANDIA_BUILD_DOCS` (bool; default: OFF): use Doxygen to build the code documentation
- `CANDIA_BUILD_MANUSCRIPTS` (bool; default: OFF): use LaTeX to build the manuscripts
- `CANDIA_BUILD_EXAMPLES` (bool; default: YES): build the two examples in the `examples/` directory
- `CANDIA_BUILD_TESTS` (bool; default: OFF): build the "tests" in the `tests/` directory. These aren't traditional tests, but rather a more extensive set of C++ files to perform e.g. benchmarking. Basic functionality is already in the examples.

Running `cmake --install <build-dir>` will perform standard installation to the chosen directory, and will also spit out a `candiaConfig.cmake` file -- see [here](#usage-from-other-cmake-projects) for what to do with it in other CMake projects.


## Usage

Once built, there will be a directory named `examples` in which there are some executables, data files, and other auxiliary files that all demonstrate a lot of the basic functionality, the source files for which are in the `examples` directory in the root of the repository.

- `evolve_dglap.cpp`: performs the evolution with the arguments passed in to the executable. Spits out a data file containing all of the resultant distributions.
- `read_table.cpp`: accepts a data file on the command line and spits out a PDF file containing a table of the results built with LaTeX. Requires `pdflatex` to be available on the command line.

Running the an executable with no arguments will indicate how to use each one. `evolve_dglap.cpp` will also contain comments on the library and what is going on.

## Usage from Other CMake Projects

As mentioned in [Compiling](#compiling), installing the project will provide a `candiaConfig.cmake` file. In a `CMakeLists.txt` file in another project, one needs only write

```cmake
enable_language(Fortran)
add_executable(main ...)
# ...
find_package(candia REQUIRED)
target_include_directories(main PUBLIC ${CANDIA_INCLUDE_DIR})
target_link_libraries(main PRIVATE ${CANDIA_LIBRARIES})
```

The call to `enable_language` is required because `Candia-v2` uses some Fortran files to do e.g. harmonic polylogarithms, and this function call sets up the required compiler/linker flags.

CMake will need to know where to find the `candiaConfig.cmake` file, which can be achieved either by appending the installation prefix to the `CMAKE_PREFIX_PATH` variable or specifying the variable `candia_DIR` to point to `<candia-prefix>/lib/cmake/candia`.

## Changelog

We will now keep a Changelog file for each release/update to the main branch, found in [Changelog.md](Changelog.md). Some minor info about what changed per release can be found on the Releases page.

## Attributions and License

This code is available under the GPLv3 license, distributed here as [LICENSE](LICENSE). We also are grateful for the code and routines provided for the NNLO and approximate N3LO splitting functions, references to which are provided in the public GitHub repository [here](https://github.com/svenolafmoch/Conformal-EIC), licensed also under the GPLv3 license. We are also grateful for the 3-loop operator matrix elements, provided via the `libome` library. The original repository is [here](https://gitlab.com/libome/libome), and a fork used in `candia-v2` is given [here](https://github.com/champso1/libome-fork), also available under the GPLv3 license, with appropriate citations required for usage given in its CITATION file.

## Contacts

Inquiries can be directed to either Casey Hampson at [champso1@students.kennesaw.edu](mailto:champso1@students.kennesaw.edu) or [casey@thehampsons.us](mailto:casey@thehampsons.us) or Dr. Marco Guzzi at [mguzzi@kennesaw.edu](mailto:mguzzi@kennesaw.edu)



