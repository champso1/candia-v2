add_executable(p3ns-exact p3ns-exact.cpp)
target_include_directories(p3ns-exact PUBLIC)
target_link_libraries(p3ns-exact PUBLIC candia ome)

file(COPY
  plot.sh
  p3nsp-approx-vs-exact-qns.gplt
  p3nsp-approx-vs-exact-lm.gplt
  p3nsm-approx-vs-exact-qns.gplt
  p3nsm-approx-vs-exact-lm.gplt
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
