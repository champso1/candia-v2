cmake_minimum_required(VERSION 3.20)
project(candia_examples LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(evolve_dglap evolve_dglap.cpp)
target_include_directories(evolve_dglap PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(evolve_dglap PUBLIC candia)

add_executable(read_table read_table.cpp)
target_include_directories(read_table PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(read_table PUBLIC candia)

# GCC places std::print functionality inside its experimental runtime library for windows
if (WIN32)
  target_link_libraries(evolve_dglap PRIVATE stdc++exp)
  target_link_libraries(read_table PRIVATE stdc++exp)
endif()
