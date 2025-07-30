cmake_minimum_required(VERSION 3.30)
project(candia_tests_cppstuff LANGUAGES CXX)

# enforce C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# make the executable, link with candia, and use the original (not installed) include dir
add_executable(cpp-test main.cpp)
