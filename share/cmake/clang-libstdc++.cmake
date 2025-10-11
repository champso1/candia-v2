set(CMAKE_C_COMPILER clang CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER clang++ CACHE STRING "" FORCE)

set(CMAKE_CXX_FLAGS_INIT "-stdlib=libstdc++")
set(CMAKE_Fortran_FLAGS_INIT "-std=legacy -w")
