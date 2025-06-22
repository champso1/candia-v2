cmake_minimum_required(VERSION 3.20)
project(candia_tests_better LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(main main.cpp)
target_compile_options(main PUBLIC -O2)

list(APPEND CMAKE_PREFIX_PATH "../../install/lib/cmake/")
find_package(Candia REQUIRED)
target_include_directories(main PUBLIC ${CANDIA_INCLUDE_DIRS})
target_link_libraries(main PUBLIC ${CANDIA_LIBRARIES})
