# study programs: they produce tables and plots for a person to look at, and are not run by ctest
# p3ns-exact and plot-generation only need the library; the others need LHAPDF
add_subdirectory(p3ns-exact)
add_subdirectory(plot-generation)

if (CANDIA_WITH_LHAPDF AND LHAPDF_FOUND)
  add_subdirectory(main)
  add_subdirectory(lhapdf-datafile)
  add_subdirectory(pdflumiunc)
else()
  message(STATUS "candia-v2: studies/main, studies/lhapdf-datafile and studies/pdflumiunc need LHAPDF (CANDIA_WITH_LHAPDF=ON); skipping them")
endif()
