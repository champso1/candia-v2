cmake_minimum_required(VERSION 3.30)
project(candia_tests_grid LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(grid grid.cpp)
set(aux
  plot.gplt
  script.py
)
file(COPY ${aux} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (WIN32)
  target_link_libraries(grid PUBLIC stdc++exp)
endif ()
