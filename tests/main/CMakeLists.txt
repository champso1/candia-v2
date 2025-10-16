cmake_minimum_required(VERSION 3.20)
project(candia_tests_main LANGUAGES CXX)

# enforce C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(evolve evolve.cpp)
target_include_directories(evolve PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(evolve PUBLIC candia)

add_executable(read read.cpp)
target_include_directories(read PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(read PUBLIC candia)

add_executable(compare compare.cpp)
target_include_directories(compare PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(compare PUBLIC candia)

add_executable(plot plot.cpp)
target_include_directories(plot PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(plot PUBLIC candia)


# must link with this stupid piece of shit library on windows
if (WIN32)
  target_link_libraries(read PUBLIC stdc++exp)
  target_link_libraries(plot PUBLIC stdc++exp)
endif()

set(resources
  VFNS-MSHT-FHMRUVV.dat
  plot.gplt
)

file(COPY ${resources} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
