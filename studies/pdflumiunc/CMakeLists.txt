add_executable(pdflumiunc pdflumiunc.cpp)
target_link_libraries(pdflumiunc PRIVATE candia LHAPDF::LHAPDF)

add_executable(evolvepdfset evolvepdfset.cpp)
target_link_libraries(evolvepdfset PRIVATE candia LHAPDF::LHAPDF)

candia_copy_to_build_dir(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/plot.gplt ${CMAKE_CURRENT_SOURCE_DIR}/infofile.in)
