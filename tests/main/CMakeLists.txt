set(all_names "evolve" "evolve-qed" "read" "compare" "ratio" "precision" "plot" "scaleratio" "test-all")

set(evolve_sources evolve.cpp)
set(evolve-qed_sources evolve-qed.cpp util.cpp)
set(read_sources read.cpp util.cpp)
set(compare_sources compare.cpp util.cpp)
set(ratio_sources ratio.cpp util.cpp)
set(precision_sources precision.cpp util.cpp)
set(plot_sources plot.cpp)
set(scaleratio_sources scaleratio.cpp util.cpp)
set(test-all_sources test-all.cpp util.cpp)

set(names "")
foreach(name ${all_names})
  candia_test_enabled(${name} want)
  if(want)
    list(APPEND names ${name})
    add_executable(${name} ${${name}_sources})
  endif()
endforeach()

foreach(name ${names})
  target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(${name} PRIVATE candia)
  target_sources(${name} PRIVATE
    FILE_SET main_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES util.hpp)
endforeach()

set(resources
  evolve.config.yaml.in
  plot-scripts
  hoppet-data
  pegasus-data
  tex-table
  tex-table2
  benchmark-data
  backup.sh
)
file(COPY ${resources} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(tests_main DEPENDS ${names})
