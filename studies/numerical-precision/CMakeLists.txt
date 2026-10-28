add_executable(numerical-precision numerical-precision.cpp)
target_link_libraries(numerical-precision PRIVATE candia)

candia_copy_to_build_dir(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/plot.gplt)
