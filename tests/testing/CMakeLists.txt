add_executable(testing testing.cpp)
target_include_directories(testing PUBLIC ${LHAPDF_INCLUDE_DIR})
target_link_libraries(testing PUBLIC candia GSL::gsl ${LHAPDF_LIBRARIES})
