add_executable(lhapdf-datafile lhapdf-datafile.cpp)
target_link_libraries(lhapdf-datafile PRIVATE candia)

candia_copy_to_build_dir(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/infofile.in ${CMAKE_CURRENT_SOURCE_DIR}/plot.gplt)
