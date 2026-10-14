add_library(studies_main_util STATIC util.cpp)
target_sources(studies_main_util PUBLIC
  FILE_SET main_headers
  TYPE HEADERS
  BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
  FILES util.hpp)
target_link_libraries(studies_main_util PUBLIC candia LHAPDF::LHAPDF)

set(studies_main_executables
  evolve
  read
  compare
  ratio
  precision
  scaleratio)
foreach(target IN LISTS studies_main_executables)
  add_executable(${target} ${target}.cpp)
  target_link_libraries(${target} PRIVATE studies_main_util)
endforeach()
add_custom_target(studies_main DEPENDS ${studies_main_executables})

file(GLOB_RECURSE resources CONFIGURE_DEPENDS plot-scripts/* tex-table/*)
candia_copy_to_build_dir(${CMAKE_CURRENT_SOURCE_DIR} ${resources} ${CMAKE_CURRENT_SOURCE_DIR}/backup.sh)
# the benchmark tables, so they can be passed as e.g. hoppet/out-nlo.dat
file(GLOB benchmark_tables CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/data/hoppet/*.dat ${PROJECT_SOURCE_DIR}/data/pegasus/*.dat)
candia_copy_to_build_dir(${PROJECT_SOURCE_DIR}/data ${benchmark_tables})
