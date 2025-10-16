cmake_minimum_required(VERSION 3.30)
project(candia_v2)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

# default the installation to the project directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)

# default to release build
# set(CMAKE_BUILD_TYPE Release)

# enforce C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# grab all source files
file(GLOB sources ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/include/*.hpp)

add_library(candia STATIC ${sources})


option(SANITIZATION "Use undefined/address sanitization (only on UNIX systems!)" off)

if (SANITIZATION)
  target_link_libraries(candia PUBLIC asan ubsan gfortran)
  target_compile_options(candia PRIVATE -Wall -Wextra -fsanitize=undefined,address -fPIC)
else()
  #target_link_libraries(candia PUBLIC gfortran)
  target_compile_options(candia PRIVATE -Wall -Wextra -fPIC)
endif()

# remove some extraneous warnings that I do not care about
target_compile_options(candia PRIVATE -Wno-unused-private-field)


# compile the fortran code for harmonic polylogs
enable_language(Fortran)
add_library(hplog STATIC src/fortran/hplog.f)
target_link_libraries(hplog PUBLIC gfortran)

# compile the fortran code for 3-loop splitting functions
enable_language(Fortran)
add_library(p2splitfuncs STATIC src/fortran/xpij2p.f src/fortran/xpns2p.f)
target_link_libraries(p2splitfuncs PUBLIC gfortran)
target_compile_options(p2splitfuncs PRIVATE -std=legacy -w)

# link these with main candia
target_link_libraries(candia PUBLIC hplog p2splitfuncs)
# on unix, we need pthread
if (UNIX)
  target_compile_options(candia PUBLIC -pthread)
  target_link_libraries(candia PUBLIC pthread)
endif()

# include directories,
target_include_directories(candia PUBLIC ${CMAKE_SOURCE_DIR}/include)


# compile tests as well
add_subdirectory(tests)


# --------------------------
#   Doxygen/Documentation
# --------------------------

# set default option to no docs
option(BUILD_DOCS "Build documentation" OFF)

if(BUILD_DOCS)
  find_package(Doxygen
    REQUIRED dot)
  if (DOXYGEN_FOUND)
	set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs)
  set(DOXYGEN_GENERATE_LATEX YES)
	set(DOXYGEN_GENERATE_TREEVIEW YES)
  set(DOXYGEN_QUIET YES)
  set(DOXYGEN_WARN_IF_UNDOCUMENTED NO)
	doxygen_add_docs(doxygen
	  ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src
	  ALL
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif(DOXYGEN_FOUND)
endif(BUILD_DOCS)

