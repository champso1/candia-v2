option(BUILD_EXAMPLES "Build the examples" YES)
if (BUILD_EXAMPLES)
  add_executable(evolve_dglap evolve_dglap.cpp)
  target_include_directories(evolve_dglap PUBLIC ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(evolve_dglap PUBLIC candia)

  add_executable(read_table read_table.cpp)
  target_include_directories(read_table PUBLIC ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(read_table PUBLIC candia)

  set(aux
	${CMAKE_CURRENT_SOURCE_DIR}/tex-table)
  file(COPY ${aux} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
