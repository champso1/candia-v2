cmake_minimum_required(VERSION 3.30)
project(candia_tests_function_grid LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(function-grid function-grid.cpp)
target_include_directories(function-grid PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(function-grid PRIVATE candia)

add_executable(function-grid-new function-grid-new.cpp)
target_include_directories(function-grid-new PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(function-grid-new PRIVATE candia)

# to use std::print on Windows
if (WIN32)
  target_link_libraries(function-grid PRIVATE stdc++exp)
  target_link_libraries(function-grid-new PRIVATE stdc++exp)
endif ()
