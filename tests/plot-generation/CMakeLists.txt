cmake_minimum_required(VERSION 3.30)
project(candia_tests_plot_generation LANGUAGES CXX)

# enforce C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# make the executable, link with candia, and use the original (not installed) include dir
add_executable(plot-generation plot-generation.cpp)
target_link_libraries(plot-generation PUBLIC candia)

file(GLOB gnuplot_files ${CMAKE_SOURCE_DIR}/tests/plot-generation/*.gplt)
file(GLOB candiav1_files ${CMAKE_SOURCE_DIR}/tests/plot-generation/candia-v1/*.dat)
file(COPY ${gnuplot_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${candiav1_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/candia-v1)

