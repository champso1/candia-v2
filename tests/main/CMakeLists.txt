add_library(tests_main_util STATIC util.cpp)
target_sources(tests_main_util PUBLIC
  FILE_SET main_headers
  TYPE HEADERS
  BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
  FILES util.hpp)
target_include_directories(tests_main_util PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(tests_main_util PUBLIC candia)

set(tests_main_executables
  evolve
  read
  compare
  ratio
  precision
  plot
  scaleratio
  test-all)
foreach(target IN LISTS tests_main_executables)
  add_executable(${target} ${target}.cpp)
  target_link_libraries(${target} PRIVATE tests_main_util)
endforeach()

set(resources
  plot-scripts
  hoppet-data
  pegasus-data
  tex-table
  benchmark-data
  backup.sh
)
file(COPY ${resources} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(tests_main DEPENDS evolve read compare ratio precision scaleratio test-all)
