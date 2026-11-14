add_executable(subtraction-pdfs subtraction-pdfs.cpp)
target_link_libraries(subtraction-pdfs PRIVATE candia LHAPDF::LHAPDF)

candia_copy_to_build_dir(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/plot.py ${CMAKE_CURRENT_SOURCE_DIR}/plot.gplt)
