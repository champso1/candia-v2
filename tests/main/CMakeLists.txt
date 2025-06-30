cmake_minimum_required(VERSION 3.20)
project(candia_tests_main LANGUAGES CXX)

# enforce C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# make the executable, link with candia, and use the original (not installed) include dir
add_executable(main main.cpp)
target_include_directories(main PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(main PUBLIC candia)
