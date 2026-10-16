add_executable(iteration-precision iteration-precision.cpp)
target_link_libraries(iteration-precision PRIVATE candia)

file(GLOB resources CONFIGURE_DEPENDS *.gplt)
candia_copy_to_build_dir(${CMAKE_CURRENT_SOURCE_DIR} ${resources})
