add_executable(pdflumiunc pdflumiunc.cpp)
target_include_directories(pdflumiunc PUBLIC ${LHAPDF_INCLUDE_DIR})
target_link_libraries(pdflumiunc PUBLIC candia ${LHAPDF_LIBRARIES})

add_executable(evolvepdfset evolvepdfset.cpp)
target_include_directories(evolvepdfset PUBLIC ${LHAPDF_INCLUDE_DIR})
target_link_libraries(evolvepdfset PUBLIC candia ${LHAPDF_LIBRARIES})

file(COPY plot.gplt infofile.in DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
