add_executable(plot-generation plot-generation.cpp)
target_link_libraries(plot-generation PRIVATE candia)

file(GLOB resources CONFIGURE_DEPENDS plot-scripts/*)
candia_copy_to_build_dir(${CMAKE_CURRENT_SOURCE_DIR} ${resources})
