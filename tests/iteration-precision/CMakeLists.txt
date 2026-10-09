add_executable(iteration-precision iteration-precision.cpp)
target_link_libraries(iteration-precision PUBLIC candia LHAPDF::LHAPDF libome::libome)

file(COPY plot.gplt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


