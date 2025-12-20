# Candia-v2

`Candia-v2` is a new and improved version of `Candia`, an algorithm/package to numerical solve the DGLAP evolution equations up to approximate next-to-next-to-next-to leading order (aN3LO) accuracy. `Candia` was originally written primarily in C and was only capable of evolution up to next-to-next-to (NNLO) accuracy. With calculations of the four-loop splitting functions and operator matrix elements being available in an approximate form, `Candia-v2` was developed in C++ in order to implement these functions and make many optimizations not present in its C counterpart.

## Building

`Candia-v2` depends on `libome`, which itself depends on `GSL`. This means that Windows support is limited. If using Windows, WSL or the MSYS2 suite is recommended. Further, at the moment, `candia-v2` depends on C++23 features, along with a fortran compiler such as `gfortran`.

Once these requirements are satisfied, building can be done with standard CMake. For ease of use, we have provided a `CMakePresets.json` which allows the user to run, for instance,

```
cmake --preset={gcc,clang}-{debug,release}-{linux,win}
```

to set up all required environment variables to build the library in debug or release mode using either `clang` or the GNU compilers on a Linux machine or Windows machine. This places build files in `bin`, which can then be compiled with `make`, `ninja`, or whatever build-tool one chooses in the configure step.
