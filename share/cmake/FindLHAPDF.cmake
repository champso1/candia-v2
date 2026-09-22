# - Find LHAPDF
#
# Finds the LHAPDF C++ library and headers.
#
# Hints/Search paths:
#   LHAPDF_ROOT          - CMake variable or environment variable pointing to the install prefix
#   CANDIA_LHAPDF_DIR    - Custom project hint
#
# Output:
#   LHAPDF_FOUND         - True if LHAPDF was found
#   LHAPDF_INCLUDE_DIRS  - Directory containing LHAPDF headers
#   LHAPDF_LIBRARIES     - Library files to link against
#   LHAPDF_VERSION       - Version string reported by lhapdf-config or pkg-config
#   LHAPDF::LHAPDF       - Imported target to link against

include(FindPackageHandleStandardArgs)

find_program(LHAPDF_CONFIG
  NAMES lhapdf-config
  HINTS
    ${LHAPDF_ROOT}
    ${CANDIA_LHAPDF_DIR}
  PATH_SUFFIXES bin)
if (LHAPDF_CONFIG)
  execute_process(
    COMMAND ${LHAPDF_CONFIG} --incdir
    OUTPUT_VARIABLE LHAPDF_INCLUDE_DIR
    ERROR_QUIET
    RESULT_VARIABLE _LHAPDF_CONFIG_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (NOT ${LHAPDF_INCDIR_RESULT} EQUAL 0)
     execute_process(
      COMMAND ${LHAPDF_CONFIG} --includedir
      OUTPUT_VARIABLE _LHAPDF_CONFIG_INCLUDE_DIR
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()

  execute_process(
    COMMAND ${LHAPDF_CONFIG} --libdir
    OUTPUT_VARIABLE _LHAPDF_CONFIG_LIBDIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
	ERROR_QUIET)

  execute_process(
    COMMAND ${LHAPDF_CONFIG} --version
    OUTPUT_VARIABLE _LHAPDF_CONFIG_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
	ERROR_QUIET)
endif()

find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
  pkg_check_modules(PC_LHAPDF LHAPDF QUIET)
endif()

find_path(LHAPDF_INCLUDE_DIR
  NAMES LHAPDF/LHAPDF.h
  HINTS
    ${_LHAPDF_CONFIG_INCLUDE_DIR}
    ${PC_LHAPDF_INCLUDE_DIRS}
    ${CANDIA_LHAPDF_DIR}
    ${LHAPDF_ROOT}
  PATH_SUFFIXES include)

find_library(LHAPDF_LIBRARY
  NAMES LHAPDF
  HINTS
    ${_LHAPDF_CONFIG_LIBDIR}
    ${PC_LHAPDF_LIBRARY_DIRS}
    ${CANDIA_LHAPDF_DIR}
    ${LHAPDF_ROOT}
  PATH_SUFFIXES lib lib64)

if(_LHAPDF_CONFIG_VERSION)
  set(LHAPDF_VERSION "${_LHAPDF_CONFIG_VERSION}")
elseif(PC_LHAPDF_VERSION)
  set(LHAPDF_VERSION "${PC_LHAPDF_VERSION}")
endif()

find_package_handle_standard_args(LHAPDF REQUIRED_VARS
  LHAPDF_INCLUDE_DIR LHAPDF_LIBRARY
  VERSION_VAR LHAPDF_VERSION)

if(LHAPDF_FOUND)
  set(LHAPDF_INCLUDE_DIRS "${LHAPDF_INCLUDE_DIRS}")
  set(LHAPDF_LIBRARIES "${LHAPDF_LIBRARIES}")

  if (NOT TARGET LHAPDF::LHAPDF)
	add_library(LHAPDF::LHAPDF UNKNOWN IMPORTED)
    set_target_properties(LHAPDF::LHAPDF PROPERTIES
      IMPORTED_LOCATION "${LHAPDF_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${LHAPDF_INCLUDE_DIR}"
    )

    # Propagate any extra flags discovered by pkg-config if applicable
    if(PC_LHAPDF_FOUND)
      set_property(TARGET LHAPDF::LHAPDF PROPERTY
        INTERFACE_COMPILE_OPTIONS "${PC_LHAPDF_CFLAGS_OTHER}")
	  set_property(TARGET LHAPDF::LHAPDF PROPERTY
		INTERFACE_LINK_LIBRARIES "${PC_LHAPDF_LDFLAGS_OTHER}")
	  message("")
    endif()
  endif()
endif()

mark_as_advanced(
  LHAPDF_CONFIG
  LHAPDF_INCLUDE_DIR
  LHAPDF_LIBRARY)
