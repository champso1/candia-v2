add_executable(plot-generation plot-generation.cpp)
target_include_directories(plot-generation PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(plot-generation PUBLIC candia)

add_executable(generate-splitfunc-plots generate-splitfunc-plots.cpp)
target_include_directories(generate-splitfunc-plots PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(generate-splitfunc-plots PUBLIC candia)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot-scripts DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

