set(names "evolve" "read" "compare" "ratio" "precision" "plot")
add_executable(compare compare.cpp util.cpp)
add_executable(evolve evolve.cpp)
add_executable(plot plot.cpp)
add_executable(precision precision.cpp util.cpp)
add_executable(ratio ratio.cpp)
add_executable(read read.cpp util.cpp)
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
  hoppet-data
  tex-table
  benchmark-data
  backup.sh
)
file(COPY ${resources} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(tests_main DEPENDS evolve read compare ratio precision)
