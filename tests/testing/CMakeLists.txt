add_executable(testing testing.cpp)
target_link_libraries(testing PUBLIC candia)

find_package(GSL REQUIRED)
target_link_libraries(testing PUBLIC GSL::gsl)

