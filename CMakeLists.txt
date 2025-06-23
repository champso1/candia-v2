cmake_minimum_required(VERSION 3.25)
project(candia_v2)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(CMAKE_BUILD_TYPE Debug)

# default the installation to the project directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)

# enforce C++20 (to avoid alignment issues with eigen) 
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# grab all source files
file(GLOB sources src/*.cpp include/*.hpp)

# make library with maximum optimization and all warnings
add_library(candia STATIC ${sources}) 
target_compile_options(candia PUBLIC -Wall -Wextra -O3 -lgfortran -fPIC)

# compile the fortran code for harmonic polylogs
enable_language(Fortran)
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -std=legacy")
set(CMAKE_Fortran_FLAGS_DEBUG "-O2 -std=legacy")
add_library(hplog STATIC src/hplog.f)
target_link_libraries(hplog PUBLIC gfortran)

target_link_libraries(candia PUBLIC hplog)

# include directories,
# differentiate while building versus the final installation
target_include_directories(candia PRIVATE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)


# --------------------------- #
# --- binary installation --- #
# --------------------------- #

# gives us CMAKE_INSTALL_<x>DIR
# that is set by CMAKE_INSTALL_PREFIX
include(GNUInstallDirs)

# installs the artifacts, and registers the target
# with an "export set"
# which essentially includes essential project info
install(TARGETS candia hplog
  EXPORT CandiaTargets
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# installs the export set as well
install(EXPORT CandiaTargets
  FILE CandiaTargets.cmake
  NAMESPACE Candia::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/candia
)

# these commands set up and install Cmake config files
# so that other projects can find this library
# it uses the Cmake.config.in file
# which these commands fill with necessary info
# it also will include the exported targets
include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_SOURCE_DIR}/share/cmake/Config.cmake.in
  ${CMAKE_BINARY_DIR}/CandiaConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/candia
)

# installs the above files
install(FILES
  ${CMAKE_BINARY_DIR}/CandiaConfig.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/candia
)


# the first install command only set some metadata up
# it only installs binaries and such
# so we have to manually install headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
