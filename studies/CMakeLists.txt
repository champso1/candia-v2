# study programs: they produce tables and plots for a person to look at, and are not run by ctest
# copies files next to the programs, keeping their paths relative to <base>;
# unlike file(COPY), configure_file copies a file again whenever it changes
function(candia_copy_to_build_dir base)
  foreach(file IN LISTS ARGN)
    file(RELATIVE_PATH relative_path ${base} ${file})
    configure_file(${file} ${CMAKE_CURRENT_BINARY_DIR}/${relative_path} COPYONLY)
  endforeach()
endfunction()

# main, pdflumiunc and subtraction-pdfs need LHAPDF; the others only need the library
add_subdirectory(p3ns-exact)
add_subdirectory(plot-generation)
add_subdirectory(lhapdf-datafile)
add_subdirectory(iteration-precision)
add_subdirectory(numerical-precision)
add_subdirectory(low-energy-evolutions)

if (CANDIA_WITH_LHAPDF AND LHAPDF_FOUND)
  add_subdirectory(main)
  add_subdirectory(pdflumiunc)
  add_subdirectory(subtraction-pdfs)
else()
  message(STATUS "candia-v2: studies/main, studies/pdflumiunc and studies/subtraction-pdfs need LHAPDF (CANDIA_WITH_LHAPDF=ON); skipping them")
endif()
