cmake_minimum_required(VERSION 3.20)
project(candia_tests_splitting_functions LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(splitting-function-test splitting-function-test.cpp)
target_include_directories(splitting-function-test PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(splitting-function-test PUBLIC candia)

set(aux_files
  ns2.gplt
  ns3.gplt
  ps2.gplt
  ps3.gplt
  ps.gplt
)
file(COPY ${aux_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
