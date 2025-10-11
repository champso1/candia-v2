cmake_minimum_required(VERSION 3.25)
project(candia_tests_gsl_convolution LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(gsl-convolution gsl-convolution.cpp)
target_include_directories(gsl-convolution PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(gsl-convolution PUBLIC candia)
