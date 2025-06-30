cmake_minimum_required(VERSION 3.25)
project(candia_v2)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

# include the toolchain file
set(CMAKE_TOOLCHAIN_FILE ${CMAKE_SOURCE_DIR}/share/cmake/toolchain.cmake)

# default the installation to the project directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)

# enforce C++20 and c99 
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# grab all source files
file(GLOB sources ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/include/*.hpp)

add_library(candia STATIC ${sources})
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-fPIC")
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall -Wextra -fPIC")
target_link_libraries(candia PUBLIC gfortran)

# compile the fortran code for harmonic polylogs
enable_language(Fortran)
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -std=legacy -w")
set(CMAKE_Fortran_FLAGS_DEBUG "-O3 -std=legacy -w")
add_library(hplog STATIC src/hplog.f)
target_link_libraries(hplog PUBLIC gfortran)

# link this with main candia
target_link_libraries(candia PUBLIC hplog)

# include directories,
# differentiate while building versus the final installation
target_include_directories(candia PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

add_subdirectory(tests)


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
