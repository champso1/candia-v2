add_subdirectory(main)

option(BUILD_ALL_TESTS "Build all tests, i.e. those other than the main" off)
if(BUILD_ALL_TESTS)
  add_subdirectory(candiav1-comparison)
  add_subdirectory(cpp-test)
  add_subdirectory(splitting-function-test)
  add_subdirectory(plot-generation)
endif()
