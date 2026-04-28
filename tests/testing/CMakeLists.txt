add_executable(testing testing.cpp)
target_sources(testing PUBLIC
  FILE_SET candia_test_testing_headers
  TYPE HEADERS
  BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
  FILES arraygrid2.hpp)
  
target_include_directories(testing PUBLIC ${LHAPDF_INCLUDE_DIR})
target_link_libraries(testing PUBLIC candia ${LHAPDF_LIBRARIES})
