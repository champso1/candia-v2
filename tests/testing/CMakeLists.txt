cmake_minimum_required(VERSION 3.30)
project(candia_tests_testing LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# make the executable, link with candia, and use the original (not installed) include dir
add_executable(testing testing.cpp)
target_link_libraries(testing PUBLIC candia)

find_package(GSL REQUIRED)
target_link_libraries(testing PUBLIC GSL::gsl GSL::gslcblas)

