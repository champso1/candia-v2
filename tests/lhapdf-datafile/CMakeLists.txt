add_executable(lhapdf-datafile lhapdf-datafile.cpp)
target_include_directories(lhapdf-datafile PUBLIC ${LHAPDF_INCLUDE_DIR})
target_link_libraries(lhapdf-datafile PUBLIC candia GSL::gsl ${LHAPDF_LIBRARIES} ome)

set(lhapdf-datafile_anc
  ${CMAKE_CURRENT_SOURCE_DIR}/run.ps1
  ${CMAKE_CURRENT_SOURCE_DIR}/infofile.in)
file(COPY ${lhapdf-datafile_anc} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
