cmake_minimum_required(VERSION 3.30)
project(candia_v2)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# build all vendor projects first before touching candia
add_subdirectory(vendor)

set(sources
  src/AlphaS.cpp
  src/Candia.cpp
  src/Distribution.cpp
  src/Evolution.cpp
  src/EvolutionThread.cpp
  src/Expression.cpp
  src/FuncArrayGrid.cpp
  src/Grid.cpp
  src/HFT.cpp
#  src/HPLog.cpp
  src/Math.cpp
  src/OperatorMatrixElements.cpp
  src/RecRel.cpp
  src/SpecialFuncs.cpp
  src/SplittingFns.cpp

  src/fortran/hplog.f
  src/fortran/xpij2p.f
  src/fortran/xpns2p.f

  include/Candia-v2/AlphaS.hpp
  include/Candia-v2/Candia.hpp
  include/Candia-v2/Common.hpp
  include/Candia-v2/Distribution.hpp
  include/Candia-v2/Expression.hpp
  include/Candia-v2/FuncArrayGrid.hpp
  include/Candia-v2/Grid.hpp
#  include/Candia-v2/HPLog.hpp
  include/Candia-v2/Math.hpp
  include/Candia-v2/OperatorMatrixElements.hpp
  include/Candia-v2/SpecialFuncs.hpp
  include/Candia-v2/SplittingFn.hpp
)

enable_language(Fortran)
add_library(hplog STATIC src/fortran/hplog.f)
target_link_libraries(hplog PUBLIC gfortran)

add_library(candia STATIC)
target_sources(candia PUBLIC ${sources})
target_link_libraries(candia PUBLIC hplog ome)
target_include_directories(candia PUBLIC ${CMAKE_SOURCE_DIR}/include)
if (UNIX)
  target_compile_options(candia PUBLIC -pthread)
  target_link_libraries(candia PUBLIC pthread)
endif()



add_subdirectory(examples)
add_subdirectory(tests)
