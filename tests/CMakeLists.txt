# regression tests, run with `ctest` from the build directory

set(CANDIA_HOPPET_DATA ${PROJECT_SOURCE_DIR}/data/hoppet)
set(CANDIA_PEGASUS_DATA ${PROJECT_SOURCE_DIR}/data/pegasus)

foreach(test IN ITEMS benchmark p3_consistency concurrency api)
  add_executable(regression_${test} ${test}.cpp)
  target_link_libraries(regression_${test} PRIVATE candia)
endforeach()

# compares against the Les Houches benchmark tables at Q = 100 GeV,
# with 13 iterations and trunc_idx 10 (the settings of the old test-all program);
# at N3LO with the approximate non-singlet P3, as used for the hoppet table
# tolerances: LO/NLO agree to ~6e-5, NNLO/N3LO currently to ~8e-4
function(candia_benchmark_test name order mur2_muf2 truncated table tolerance)
  add_test(NAME benchmark.${name}
    COMMAND regression_benchmark ${order} 13 10 ${mur2_muf2} ${truncated} ${table} ${tolerance})
endfunction()

candia_benchmark_test(lo.exact         0 1.0 0 ${CANDIA_HOPPET_DATA}/out-lo.dat        1e-4)
candia_benchmark_test(nlo.exact        1 1.0 0 ${CANDIA_HOPPET_DATA}/out-nlo.dat       1e-4)
candia_benchmark_test(nlo.truncated    1 1.0 1 ${CANDIA_HOPPET_DATA}/out-nlo.dat       1e-4)
candia_benchmark_test(nnlo.exact       2 1.0 0 ${CANDIA_HOPPET_DATA}/out-nnlo.dat      1e-3)
candia_benchmark_test(nnlo.truncated   2 1.0 1 ${CANDIA_HOPPET_DATA}/out-nnlo.dat      1e-3)
candia_benchmark_test(n3lo.truncated   3 1.0 1 ${CANDIA_HOPPET_DATA}/out-n3lo.dat      1e-3)
candia_benchmark_test(n3lo.exact       3 1.0 0 ${CANDIA_HOPPET_DATA}/out-n3lo.dat      1e-3)
candia_benchmark_test(lo.mu20          0 2.0 0 ${CANDIA_PEGASUS_DATA}/out-lo-mu20.dat  1e-4)
candia_benchmark_test(lo.mu05          0 0.5 0 ${CANDIA_PEGASUS_DATA}/out-lo-mu05.dat  1e-4)
candia_benchmark_test(nlo.mu20         1 2.0 0 ${CANDIA_PEGASUS_DATA}/out-nlo-mu20.dat 1e-4)
candia_benchmark_test(nlo.mu05         1 0.5 0 ${CANDIA_PEGASUS_DATA}/out-nlo-mu05.dat 1e-4)

add_test(NAME p3_consistency COMMAND regression_p3_consistency)
add_test(NAME concurrency COMMAND regression_concurrency)
add_test(NAME api COMMAND regression_api)

if (CANDIA_WITH_LHAPDF AND LHAPDF_FOUND)
  add_executable(regression_lhapdf_roundtrip lhapdf_roundtrip.cpp)
  target_link_libraries(regression_lhapdf_roundtrip PRIVATE candia LHAPDF::LHAPDF)
  set(roundtrip_dir ${CMAKE_CURRENT_BINARY_DIR}/lhapdf_roundtrip)
  file(MAKE_DIRECTORY ${roundtrip_dir})
  add_test(NAME lhapdf_roundtrip
    COMMAND regression_lhapdf_roundtrip ${PROJECT_SOURCE_DIR}/examples/infofile.in
    WORKING_DIRECTORY ${roundtrip_dir})
endif()

# every evolution test uses all cores; running them in parallel only oversubscribes
get_property(candia_regression_tests DIRECTORY PROPERTY TESTS)
set_tests_properties(${candia_regression_tests} PROPERTIES RUN_SERIAL TRUE)
