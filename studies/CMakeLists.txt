# study programs: they produce tables and plots for a person to look at, and are not run by ctest
# copies files next to the programs, keeping their paths relative to <base>;
# unlike file(COPY), configure_file copies a file again whenever it changes
function(candia_copy_to_build_dir base)
  foreach(file IN LISTS ARGN)
    file(RELATIVE_PATH relative_path ${base} ${file})
    configure_file(${file} ${CMAKE_CURRENT_BINARY_DIR}/${relative_path} COPYONLY)
  endforeach()
endfunction()

# main and pdflumiunc need LHAPDF; the others only need the library
add_subdirectory(p3ns-exact)
add_subdirectory(plot-generation)
add_subdirectory(lhapdf-datafile)

if (CANDIA_WITH_LHAPDF AND LHAPDF_FOUND)
  add_subdirectory(main)
  add_subdirectory(pdflumiunc)
else()
  message(STATUS "candia-v2: studies/main and studies/pdflumiunc need LHAPDF (CANDIA_WITH_LHAPDF=ON); skipping them")
endif()
