add_executable(plot-generation plot-generation.cpp)
target_include_directories(plot-generation PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(plot-generation PUBLIC candia)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot-scripts DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

