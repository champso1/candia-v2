option(CANDIA_BUILD_TESTS "Build the 'tests', which are effectively just a more exhaustive set of executables for testing the evolution." OFF)

set(CANDIA_TEST_TARGETS "" CACHE STRING
  "Semicolon-separated list of test executables to build, e.g. -DCANDIA_TEST_TARGETS=\"evolve;read\". Empty (default) builds all available tests.")

# Sets ${out_var} to TRUE if 'name' should be built: either CANDIA_TEST_TARGETS
# is empty (build everything) or it explicitly lists 'name'.
function(candia_test_enabled name out_var)
  if(NOT CANDIA_TEST_TARGETS OR name IN_LIST CANDIA_TEST_TARGETS)
    set(${out_var} TRUE PARENT_SCOPE)
  else()
    set(${out_var} FALSE PARENT_SCOPE)
  endif()
endfunction()

if (CANDIA_BUILD_TESTS)
  add_subdirectory(main)
  add_subdirectory(plot-generation)
  if (CANDIA_WITH_LHAPDF)
	add_subdirectory(testing)
	add_subdirectory(p3ns-exact)
	add_subdirectory(lhapdf-datafile)
	add_subdirectory(pdflumiunc)
  endif()
endif()
