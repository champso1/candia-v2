cmake_minimum_required(VERSION 3.20)
project(candia_v2 VERSION 1.5.0)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message("Defaulting to release build.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
	"Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(BUILD_MANUSCRIPTS "Use pdflatex/bibtex to compile the available manuscripts" OFF)
if (BUILD_MANUSCRIPTS)
  message("candia-v2: building manuscripts")
  add_subdirectory(manuscripts)
endif()

# build all vendor projects first before touching candia
add_subdirectory(vendor)

set(sources
  src/AlphaS.cpp
  src/ArrayGrid.cpp
  src/Candia.cpp
  src/Distribution.cpp
#  src/Evolution.cpp
#  src/EvolutionThread.cpp
  src/Expression.cpp
  src/Grid.cpp
  src/HFT.cpp
  src/LHAPDFDataFile.cpp
#  src/LHAPDFDistribution.cpp
  src/Math.cpp
  src/OperatorMatrixElements.cpp
  src/RecRel.cpp
  src/SpecialFuncs.cpp
  src/SplittingFns.cpp)

set(headers
  include/Candia-v2/AlphaS.hpp
  include/Candia-v2/ArrayGrid.hpp
  include/Candia-v2/Candia.hpp
  include/Candia-v2/Common.hpp
  include/Candia-v2/Distribution.hpp
  include/Candia-v2/Expression.hpp
  include/Candia-v2/Grid.hpp
  include/Candia-v2/LHAPDFDataFile.hpp
#  include/Candia-v2/LHAPDFDistribution.hpp
  include/Candia-v2/Math.hpp
  include/Candia-v2/OperatorMatrixElements.hpp
  include/Candia-v2/Options.hpp
  include/Candia-v2/SpecialFuncs.hpp
  include/Candia-v2/SplittingFn.hpp)

enable_language(Fortran)
add_library(hplog STATIC src/fortran/hplog5.f)

set(nnlo_splitfuncs_fortran_src
  src/fortran/mvv_p2/xpij2p.f
  src/fortran/mvv_p2/xpns2p.f)
add_library(nnlo_splitfuncs_fortran STATIC ${nnlo_splitfuncs_fortran_src})

set(n3lo_splitfuncs_fortran_src
  src/fortran/mvv_p3/xpgg3a.f90
  src/fortran/mvv_p3/xpgg3a_2410.f90
  src/fortran/mvv_p3/xpgq3a.f90
  src/fortran/mvv_p3/xpgq3a_2404.f90
  src/fortran/mvv_p3/xpgq3a_2512.f90
  src/fortran/mvv_p3/xpns3m.f90
  src/fortran/mvv_p3/xpns3p.f90
  src/fortran/mvv_p3/xpns3s.f90
  src/fortran/mvv_p3/xpps3a.f90
  src/fortran/mvv_p3/xpqg3a.f90)
add_library(n3lo_splitfuncs_fortran STATIC ${n3lo_splitfuncs_fortran_src})

find_package(GSL REQUIRED)

add_library(candia STATIC)
target_sources(candia PRIVATE ${sources})
target_sources(candia PUBLIC
  FILE_SET cxx_headers
  TYPE HEADERS
  BASE_DIRS ${CMAKE_SOURCE_DIR}/include
  FILES ${headers})

option(WITH_LHAPDF "Compile in LHAPDF interface" OFF)
option(LHAPDF_ROOT_DIR "Path to the installation root/prefix of LHAPDF." "/usr")
if (WITH_LHAPDF)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/share/cmake")
  find_package(LHAPDF REQUIRED)
  message("candia-v2: Compiling with LHAPDF support.")
  target_include_directories(candia PUBLIC ${LHAPDF_INCLUDE_DIR})
  target_link_libraries(candia PUBLIC ${LHAPDF_LIBRARIES})

  target_sources(candia PUBLIC
    FILE_SET cxx_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_SOURCE_DIR}/include
    FILES include/Candia-v2/LHAPDFDistribution.hpp)
  target_sources(candia PRIVATE src/LHAPDFDistribution.cpp)
endif()

option(ENABLE_THREADING "Use Multi-Threading capabilities to perform the evolution." OFF)
if(ENABLE_THREADING)
  message("candia-v2: Compiling with threading support.")
  find_package(Threads REQUIRED)
  target_sources(candia PRIVATE src/EvolutionThread.cpp)
  target_link_libraries(candia PRIVATE Threads::Threads)
  target_compile_definitions(candia PRIVATE ENABLE_THREADING=1)
else()
  target_sources(candia PRIVATE src/Evolution.cpp)
  target_compile_definitions(candia PRIVATE ENABLE_THREADING=0)
endif()


target_link_libraries(candia PRIVATE hplog nnlo_splitfuncs_fortran n3lo_splitfuncs_fortran ome GSL::gsl)
target_include_directories(candia PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_include_directories(candia PUBLIC ${GSL_INCLUDE_DIRS})
add_subdirectory(examples)
add_subdirectory(tests)

# ====================
# DOXYGEN
# ====================

option(BUILD_DOCS "Build Doxygen documentation" OFF)

if(BUILD_DOCS)
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
	set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs)
	set(DOXYGEN_GENERATE_TREEVIEW YES)
	set(DOXYGEN_HAVE_DOT YES)
	doxygen_add_docs(doxygen
	  ${headers} ${sources}
	  ALL)
  endif(DOXYGEN_FOUND)
endif(BUILD_DOCS)

# ====================
# INSTALLATION
# ====================
include(GNUInstallDirs)
install(TARGETS candia hplog nnlo_splitfuncs_fortran n3lo_splitfuncs_fortran
  EXPORT candiaTargets
  FILE_SET cxx_headers DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT candiaTargets
  NAMESPACE candia::
  FILE candiaTargets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/candia)

set(INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR})
include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/share/cmake/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/candiaConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/candia
  PATH_VARS INSTALL_INCLUDE_DIR)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/candiaConfig.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/candia)

include(CPack)
