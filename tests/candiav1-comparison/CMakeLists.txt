cmake_minimum_required(VERSION 3.20)
project(candiav1_comparison)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)




add_executable(candiav2 candiav2/candiav2.cpp)
target_compile_options(candiav2 PUBLIC -Wall -Wextra)
target_include_directories(candiav2 PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(candiav2 PUBLIC hplog candia)



enable_language(Fortran)
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -std=legacy -fPIE -w")
set(CMAKE_Fortran_FLAGS_DEBUG "-O2 -std=legacy -fPIE -w")
add_library(candiav1-fortran candiav1/hplog.f
  candiav1/xpns2p.f
  candiav1/xpij2p.f
  candiav1/partonww.f
  candiav1/a02m.f)
target_link_libraries(candiav1-fortran PUBLIC gfortran)

add_executable(candiav1 candiav1/candia.c)
target_link_libraries(candiav1 PUBLIC candiav1-fortran)
target_compile_options(candiav1 PUBLIC -w)



add_executable(comparison comparison.cpp)
target_compile_options(comparison PUBLIC -Wall -Wextra)

