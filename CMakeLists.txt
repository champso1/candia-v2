cmake_minimum_required(VERSION 3.30)
project(candia_v2)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(ENABLE_THREADING "Enable multi-threading support" OFF)

if (ENABLE_THREADING)
  message("candia-v2: Compiling with threading support.")
endif()

# build all vendor projects first before touching candia
add_subdirectory(vendor)

set(sources
  src/AlphaS.cpp
  src/Candia.cpp
  src/Distribution.cpp
  src/Evolution.cpp
  src/Expression.cpp
  src/ArrayGrid.cpp
  src/Grid.cpp
  src/HFT.cpp
  src/Math.cpp
  src/OperatorMatrixElements.cpp
  src/RecRel.cpp
  src/SpecialFuncs.cpp
  src/SplittingFns.cpp

  src/fortran/hplog.f
  src/fortran/xpij2p.f
  src/fortran/xpns2p.f)
set(headers
  include/Candia-v2/AlphaS.hpp
  include/Candia-v2/Candia.hpp
  include/Candia-v2/Common.hpp
  include/Candia-v2/Distribution.hpp
  include/Candia-v2/Expression.hpp
  include/Candia-v2/ArrayGrid.hpp
  include/Candia-v2/Grid.hpp
  include/Candia-v2/Math.hpp
  include/Candia-v2/OperatorMatrixElements.hpp
  include/Candia-v2/SpecialFuncs.hpp
  include/Candia-v2/SplittingFn.hpp)

enable_language(Fortran)
add_library(hplog STATIC src/fortran/hplog.f)
target_link_libraries(hplog PUBLIC gfortran)

find_package(GSL REQUIRED)

option(LHAPDF_ROOT_DIR "Path to the installation root/prefix of LHAPDF." "/usr")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/share/cmake")
find_package(LHAPDF REQUIRED)

add_library(candia STATIC)
target_sources(candia PRIVATE ${sources})
target_sources(candia PUBLIC
  FILE_SET cxx_headers
  TYPE HEADERS
  BASE_DIRS ${CMAKE_SOURCE_DIR}/include
  FILES ${headers})

if(ENABLE_THREADING)
  find_package(Threads REQUIRED)
  target_sources(candia PRIVATE src/EvolutionThread.cpp)
  target_link_libraries(candia PRIVATE Threads::Threads)
  target_compile_definitions(candia PRIVATE ENABLE_THREADING=1)
else()
  target_compile_definitions(candia PRIVATE ENABLE_THREADING=0)
endif()

target_link_libraries(candia PRIVATE hplog ome ${LHAPDF_LIBRARIES} GSL::gsl)
target_include_directories(candia PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_include_directories(candia PUBLIC ${LHAPDF_INCLUDE_DIR})
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
install(TARGETS candia hplog
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

