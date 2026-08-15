cmake_minimum_required(VERSION 3.20)
project(candia_v2 VERSION 1.8.0)
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

if (NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)
endif()

# build all vendor projects first before touching candia
# add_subdirectory(vendor)


set(sources
  src/AlphaS.cpp
  src/Candia.cpp
  src/Distribution.cpp
  src/Evolution.cpp
  src/Expression.cpp
  src/Grid.cpp
  src/HFT.cpp
#  src/LHAPDFDistribution.cpp
  src/LHAPDFGrid.cpp
  src/LibomeInterface.cpp
  src/Math.cpp
  src/OperatorMatrixElements.cpp
  src/RecRel.cpp
  src/SpecialFuncs.cpp
  src/SplittingFns.cpp
  src/SplittingFnsQED.cpp
  src/SplittingFns_P3nsExact.cpp)

set(headers
  include/Candia-v2/AlphaS.hpp
  include/Candia-v2/ArrayGrid.hpp
  include/Candia-v2/Candia.hpp
  include/Candia-v2/Common.hpp
  include/Candia-v2/Couplings.hpp
  include/Candia-v2/Distribution.hpp
  include/Candia-v2/Expression.hpp
  include/Candia-v2/Grid.hpp
#  include/Candia-v2/LHAPDFDistribution.hpp
  include/Candia-v2/LHAPDFGrid.hpp
  include/Candia-v2/LibomeInterface.hpp
  include/Candia-v2/Math.hpp
  include/Candia-v2/OperatorMatrixElements.hpp
  include/Candia-v2/Options.hpp
  include/Candia-v2/SpecialFuncs.hpp
  include/Candia-v2/SplittingFn.hpp
  include/Candia-v2/SplittingFnQED.hpp)

enable_language(Fortran)

set(hplog_src src/fortran/hplog5.f)
set(nnlo_splitfuncs_fortran_src
  src/fortran/mvv_p2/xpij2p.f
  src/fortran/mvv_p2/xpns2p.f)
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

add_library(candia SHARED)

option(CANDIA_ENABLE_WARNINGS "Enable compilation with max warnings e.g. -Wall -Wextra and -Wpedantic")
if (CANDIA_ENABLE_WARNINGS)
  target_compile_options(candia PUBLIC -Wall -Wextra -Wpedantic)
endif()
target_sources(candia PRIVATE ${sources})
target_sources(candia PRIVATE
  ${hplog_src}
  ${nnlo_splitfuncs_fortran_src}
  ${n3lo_splitfuncs_fortran_src}
)
target_sources(candia PUBLIC
  FILE_SET cxx_headers
  TYPE HEADERS
  BASE_DIRS ${CMAKE_SOURCE_DIR}/include
  FILES ${headers})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/share/cmake")
find_package(TBB REQUIRED)
find_package(Threads REQUIRED)
find_package(libome REQUIRED)

target_link_libraries(candia PRIVATE
  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}
  Threads::Threads
  TBB::tbb
  libome::libome)
target_include_directories(candia PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

option(CANDIA_WITH_LHAPDF "Compile in LHAPDF interface" OFF)
option(CANDIA_LHAPDF_DIR "Path to the installation root/prefix of LHAPDF." "/usr/local")
if (CANDIA_WITH_LHAPDF)
  find_package(LHAPDF REQUIRED)
  message("candia-v2: Compiling with LHAPDF support.")
  target_link_libraries(candia PUBLIC LHAPDF::LHAPDF)

  # above, PUBLIC was necessary because we are "exposing" LHAPDF headers
  # and such in the public interface in LHAPDFDistribution.hpp
  # so down the line, in the installation interface,
  # users linking with candia::candia should have LHAPDF
  # found (handled by the installed FindLHAPDF.cmake script)
  # and linked automatically

  target_sources(candia PUBLIC
    FILE_SET cxx_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_SOURCE_DIR}/include
    FILES include/Candia-v2/LHAPDFDistribution.hpp)
  target_sources(candia PRIVATE src/LHAPDFDistribution.cpp)
endif()


add_subdirectory(examples)
add_subdirectory(tests)



# =========================
# NON-CODE-RELATED
# =========================
option(CANDIA_BUILD_DOCS "Build Doxygen documentation" OFF)
if(CANDIA_BUILD_DOCS)
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
	set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs)
	set(DOXYGEN_GENERATE_TREEVIEW YES)
	set(DOXYGEN_HAVE_DOT YES)
	doxygen_add_docs(doxygen
	  ${headers} ${sources}
	  ALL)
  endif()
endif()


# ====================
# INSTALLATION
# ====================
include(GNUInstallDirs)
install(TARGETS candia
  EXPORT candiaTargets
  FILE_SET cxx_headers DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT candiaTargets
  NAMESPACE candia::
  FILE candiaTargets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/candia)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/share/cmake/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/candiaConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/candia)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/candiaConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/candiaConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/candiaConfigVersion.cmake
    ${CMAKE_SOURCE_DIR}/share/cmake/FindLHAPDF.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/candia)


# ====================
# SOURCE TARBALL (CPack)
# ====================
# make better looking tarballs for github releases
set(CPACK_PACKAGE_NAME "candia-v2")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "candia-v2-${PROJECT_VERSION}-src")
set(CPACK_SOURCE_IGNORE_FILES
  "/\\\\.git$"
  "/\\\\.git/"
  "/\\\\.cache/"
  "/\\\\.vscode/"
  "/build[^/]*/"
  "/install/"
  "/bin/"
  "/manuscripts/temp/"
  "/manuscripts/cpc-src/"
  "\\\\.tar\\\\.gz$")
include(CPack)
