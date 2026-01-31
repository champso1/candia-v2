include(UseLATEX.cmake)
find_package(LATEX COMPONENTS PDFLATEX BIBER)

add_latex_document(candia-v2-n3lo-cpc.tex
    TARGET_NAME cpc)

add_latex_document(DGLAP-N3LO.tex
    BIBFILES sources.bib
    IMAGE_DIRS res/
    TARGET_NAME prd)

  set_property(TARGET prd_pdf PROPERTY COMPILER_FLAGS "-shell-escape")
  set_property(TARGET cpc_pdf PROPERTY COMPILER_FLAGS "-shell-escape")
  
file(COPY res DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(manuscripts_pdf DEPENDS cpc_pdf prd_pdf)
