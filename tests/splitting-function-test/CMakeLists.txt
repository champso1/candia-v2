cmake_minimum_required(VERSION 3.30)
project(candia_tests_splitting_functions LANGUAGES CXX)

# enforce C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# make the executable, link with candia, and use the original (not installed) include dir
add_executable(splitting-function-test splitting-function-test.cpp)
target_link_libraries(splitting-function-test PUBLIC candia)

file(GLOB gplt_files ${CMAKE_CURRENT_SOURCE_DIR}/*.gplt)
file(COPY ${gplt_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
