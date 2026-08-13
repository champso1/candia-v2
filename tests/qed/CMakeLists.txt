add_executable(qed qed.cpp)
target_link_libraries(qed PUBLIC candia)

file(COPY plot.gplt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
