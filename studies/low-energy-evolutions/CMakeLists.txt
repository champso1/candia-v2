add_executable(low-energy-evolutions low-energy-evolutions.cpp)
target_link_libraries(low-energy-evolutions PRIVATE candia)

candia_copy_to_build_dir(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/plot.gplt)
