add_executable(p3ns-exact p3ns-exact.cpp)
target_link_libraries(p3ns-exact PRIVATE candia)

file(GLOB resources CONFIGURE_DEPENDS *.gplt plot.sh)
candia_copy_to_build_dir(${CMAKE_CURRENT_SOURCE_DIR} ${resources})
