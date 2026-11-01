# regression tests, run with `ctest` from the build directory

set(CANDIA_HOPPET_DATA ${PROJECT_SOURCE_DIR}/data/hoppet)
set(CANDIA_PEGASUS_DATA ${PROJECT_SOURCE_DIR}/data/pegasus)

foreach(test IN ITEMS benchmark p3_consistency concurrency api)
  add_executable(regression_${test} ${test}.cpp)
  target_link_libraries(regression_${test} PRIVATE candia)
endforeach()

# compares against the Les Houches benchmark tables at Q = 100 GeV (data/hoppet: hoppet-v2, see data/README.md;
# data/pegasus), with 13 iterations and trunc_idx 10 (the settings of the old test-all program)
# agreement: ~1e-5 with hoppet up to NNLO, ~6e-5 with pegasus for mu_R != mu_F
function(candia_benchmark_test name order mur2_muf2 truncated table tolerance)
  add_test(NAME benchmark.${name}
    COMMAND regression_benchmark ${order} 13 10 ${mur2_muf2} ${truncated} ${table} ${tolerance})
endfunction()

candia_benchmark_test(lo.exact         0 1.0 0 ${CANDIA_HOPPET_DATA}/out-lo.dat        5e-5)
candia_benchmark_test(nlo.exact        1 1.0 0 ${CANDIA_HOPPET_DATA}/out-nlo.dat       5e-5)
candia_benchmark_test(nlo.truncated    1 1.0 1 ${CANDIA_HOPPET_DATA}/out-nlo.dat       5e-5)
candia_benchmark_test(nnlo.exact       2 1.0 0 ${CANDIA_HOPPET_DATA}/out-nnlo.dat      5e-5)
candia_benchmark_test(nnlo.truncated   2 1.0 1 ${CANDIA_HOPPET_DATA}/out-nnlo.dat      5e-5)
candia_benchmark_test(lo.mu20          0 2.0 0 ${CANDIA_PEGASUS_DATA}/out-lo-mu20.dat  1e-4)
candia_benchmark_test(lo.mu05          0 0.5 0 ${CANDIA_PEGASUS_DATA}/out-lo-mu05.dat  1e-4)
candia_benchmark_test(nlo.mu20         1 2.0 0 ${CANDIA_PEGASUS_DATA}/out-nlo-mu20.dat 1e-4)
candia_benchmark_test(nlo.mu05         1 0.5 0 ${CANDIA_PEGASUS_DATA}/out-nlo-mu05.dat 1e-4)

# N3LO against hoppet, once with the approximate non-singlet P3 (*-approx tables) and once with the exact one
# (*-exact tables). The two differ by only 4e-5 to 1.5e-4 at Q = 100 GeV, so these run converged (20 iterations,
# trunc_idx 15; with the defaults the truncation error is as large) and the tolerances sit between the agreement
# with the right table (<= 4.7e-5, 2.3e-5 and 6.5e-6) and with the other variant's (>= 1.0e-4, 1.4e-4 and 4.3e-5)
function(candia_benchmark_n3lo_test name truncated p3ns ffns_nf mur2_muf2 table tolerance)
  add_test(NAME benchmark.${name}
    COMMAND regression_benchmark 3 20 15 ${mur2_muf2} ${truncated} ${table} ${tolerance} ${ffns_nf} ${p3ns})
endfunction()

foreach(p3ns IN ITEMS approx exact)
  candia_benchmark_n3lo_test(n3lo.${p3ns}-p3ns.exact      0 ${p3ns} 0 1.0 ${CANDIA_HOPPET_DATA}/out-n3lo-${p3ns}.dat 7e-5)
  candia_benchmark_n3lo_test(n3lo.${p3ns}-p3ns.truncated  1 ${p3ns} 0 1.0 ${CANDIA_HOPPET_DATA}/out-n3lo-${p3ns}.dat 7e-5)
  candia_benchmark_n3lo_test(n3lo.${p3ns}-p3ns.ffn4.mu05  0 ${p3ns} 4 0.5 ${CANDIA_HOPPET_DATA}/out-n3lo-${p3ns}-ffn4-mu05.dat 5e-5)
  candia_benchmark_n3lo_test(n3lo.${p3ns}-p3ns.ffn4.mu20  0 ${p3ns} 4 2.0 ${CANDIA_HOPPET_DATA}/out-n3lo-${p3ns}-ffn4-mu20.dat 2e-5)
endforeach()

# fixed flavour number scheme (n_f = 4) with mu_R != mu_F, against hoppet;
# the single evolution from Q0 to 100 GeV needs 20 iterations to converge at large x; they agree to ~5e-5
function(candia_benchmark_ffns_test name order mur2_muf2 table)
  add_test(NAME benchmark.${name}
    COMMAND regression_benchmark ${order} 20 10 ${mur2_muf2} 0 ${table} 1e-4 4)
endfunction()

candia_benchmark_ffns_test(lo.ffn4.mu20   0 2.0 ${CANDIA_HOPPET_DATA}/out-lo-ffn4-mu20.dat)
candia_benchmark_ffns_test(lo.ffn4.mu05   0 0.5 ${CANDIA_HOPPET_DATA}/out-lo-ffn4-mu05.dat)
candia_benchmark_ffns_test(nlo.ffn4.mu20  1 2.0 ${CANDIA_HOPPET_DATA}/out-nlo-ffn4-mu20.dat)
candia_benchmark_ffns_test(nlo.ffn4.mu05  1 0.5 ${CANDIA_HOPPET_DATA}/out-nlo-ffn4-mu05.dat)
candia_benchmark_ffns_test(nnlo.ffn4.mu20 2 2.0 ${CANDIA_HOPPET_DATA}/out-nnlo-ffn4-mu20.dat)
candia_benchmark_ffns_test(nnlo.ffn4.mu05 2 0.5 ${CANDIA_HOPPET_DATA}/out-nnlo-ffn4-mu05.dat)

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
