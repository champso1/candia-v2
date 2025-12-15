cmake_minimum_required(VERSION 3.30)
project(candia_tests_plot_generation LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(plot-generation plot-generation.cpp)
target_link_libraries(plot-generation PUBLIC candia)

set(gnuplot-scripts
  ${CMAKE_CURRENT_SOURCE_DIR}/ratio.gplt)

file(COPY ${gnuplot-scripts} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scripts)

