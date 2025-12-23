# Candia-v2

**NOTE: `Candia-v2` is still a work-in-progress. The code works correctly and basic functionality is present, but it is not perfectly optimized, nor is the code very clean or feature-rich. We are actively working hard to expand and improve the code.**

`Candia-v2` is a new and improved version of `Candia` (cfr. A. Cafarella, M. Guzzi, C. Coriano', Comp.Phys.Comm. 179 2008, A. Cafarella, M. Guzzi, C. Coriano', Nucl.Phys.B748 2006), a computer code to numerically solve DGLAP evolution for collinear PDFs in the x-space up to next-to-next-to-next-to leading order (N^3LO) accuracy in perturbative QCD. `Candia-v2` currently uses an approximate version of the 4-loop splitting functions as the calculation of their exact analytical form is still in progress. `Candia` was originally written in C and was only capable of evolution up to next-to-next-to (NNLO) accuracy in QCD. New information on the splitting functions and operator matrix elements as well as advancements to C++ led to the development of `Candia-v2`.

## Building

### Prerequisites

`Candia-v2` depends on `libome`, which is in this repository as a submodule. When cloning this repository, either pass the option `--recurse-submodules` to `git clone`, or after the repository is cloned, run `git submodule init` followed by `git submodule update`. `libome` itself depends on GSL, which means that `candia-v2` can only be built on UNIX-like systems. Windows is supported either via WSL, the MSYS2 suite, or Cygwin (the former two are tested, but Cygwin is not). Lastly, at the moment, there are several C++23 features used throughout, such as `std::print` and the `ranges` suite. We are planning on reverting some of these things in order to have C++17 support, but this is not the case at the moment.

### Compiling

Compiling follows the standard CMake procedure:

```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

Specifying a release build will of course make the code run significantly faster.

## Usage

Once `candia-v2` is built, there will be a directory named `examples` in which there are some executables, data files, and other auxiliary files that all demonstrate a lot of the basic functionality, the source files for which are in the `examples` directory in the root of the repository.

- `evolve_dglap.cpp`: performs the evolution with the arguments passed in to the executable. Spits out a data file containing all of the resultant distributions.
- `read_table.cpp`: accepts a data file on the command line and spits out a PDF file containing a table of the results built with LaTeX. Requires `pdflatex` to be available on the command line.

Running the associated executables with no options will indicate how to use each one. `evolve_dglap.cpp` also contains comments on the code related to `candia-v2` (`read_table.cpp` contains just simple C++ code for file I/O and other non-`candia-v2` related stuff).

## Attributions and License

This code is available under the GPLv3 license, distributed here as [LICENSE](LICENSE). We also are grateful for the code and routines provided for the NNLO and approximate N^3LO splitting functions, references to which are provided in the public GitHub repository [here](https://github.com/svenolafmoch/Conformal-EIC), licensed also under the GPLv3 license. We are also grateful for the 3-loop operator matrix elements, provided via the `libome` library. The original repository is [here](https://gitlab.com/libome/libome), and a fork used in `candia-v2` is given [here](https://github.com/champso1/libome-fork), also available under the GPLv3 license, with appropriate citations required for usage given in its CITATION file.



