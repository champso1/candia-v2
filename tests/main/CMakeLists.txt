cmake_minimum_required(VERSION 3.20)
project(candia_tests_main LANGUAGES CXX)

# enforce C++ standard
set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(evolve evolve.cpp)
target_include_directories(evolve PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(evolve PUBLIC candia)

add_executable(read read.cpp)
target_include_directories(read PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(read PUBLIC candia)

add_executable(compare-hoppet compare-hoppet.cpp)
target_include_directories(compare-hoppet PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(compare-hoppet PUBLIC candia)

add_executable(ratio ratio.cpp)
target_include_directories(ratio PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(ratio PUBLIC candia)

# must link with this stupid piece of shit library on windows
# to use ordinary C++20 println features man come on
if (WIN32)
  target_link_libraries(evolve PRIVATE stdc++exp)
  target_link_libraries(read PRIVATE stdc++exp)
  target_link_libraries(compare-hoppet PRIVATE stdc++exp)
  target_link_libraries(ratio PRIVATE stdc++exp)
endif()

set(resources
  out.dat
  out-ffns.dat
  out-nnlo.dat
)

file(COPY ${resources} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
