add_executable(lhapdf-datafile lhapdf-datafile.cpp)
target_include_directories(lhapdf-datafile PUBLIC ${LHAPDF_INCLUDE_DIR})
target_link_libraries(lhapdf-datafile PUBLIC candia ${LHAPDF_LIBRARIES} ome)

add_executable(subpdfs subpdfs.cpp)
target_include_directories(subpdfs PUBLIC ${LHAPDF_INCLUDE_DIR})
target_link_libraries(subpdfs PUBLIC candia ${LHAPDF_LIBRARIES} ome)

set(lhapdf-datafile_anc
  infofile.in
  plot.gplt)
file(COPY ${lhapdf-datafile_anc} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
