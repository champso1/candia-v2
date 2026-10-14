add_executable(pdflumiunc pdflumiunc.cpp)
target_link_libraries(pdflumiunc PUBLIC candia LHAPDF::LHAPDF)

add_executable(evolvepdfset evolvepdfset.cpp)
target_link_libraries(evolvepdfset PUBLIC candia LHAPDF::LHAPDF)

file(COPY plot.gplt infofile.in DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
