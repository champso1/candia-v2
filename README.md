# Candia-v2

`Candia-v2` is a new and improved version of `Candia`, an algorithm/package to numerical solve the DGLAP evolution equations up to approximate next-to-next-to-next-to leading order (aN3LO) accuracy. `Candia` was originally written primarily in C and was only capable of evolution up to next-to-next-to (NNLO) accuracy. With calculations of the four-loop splitting functions and operator matrix elements being available in an approximate form, `Candia-v2` was developed in C++ in order to implement these functions and make many optimizations not present in its C counterpart.

## Building

`Candia-v2` depends on `libome`, which itself depends on `GSL`. This means that Windows support is limited. If using Windows, WSL is preferred. There is also the option of MSYS2, a set of tools providing a UNIX-like interface as well as many programs typically available on UNIX-like machines (such as GSL).  Additionally, `libome` requires the GNU compilers, namely `g++`, as CLang fails to compile correctly. Lastly, there are many modern C++ features used in `Candia-v2`, and thus it requires C++23 support.

Once these requirements are satisfied, building can be done with standard CMake. For ease of use, we have provided a `CMakePresets.json` which allows the user to run, for instance,

```
cmake --preset=gcc-release-linux
```

to set up all required environment variables to build the library in release mode. This places build files in `bin/gcc-release-linux`. The default build tool is `Ninja`. One then must

```
cd bin/gcc-release-linux
ninja
```

to build the library (and, by default, the examples).

## Examples (TODO)

By default, the above build procedure builds a few examples which illustrate the usage of the library. There are also a few other scripts that are compiled to provide some plotting functionality via the `gnuplot` tool as well as table generation with LaTeX.