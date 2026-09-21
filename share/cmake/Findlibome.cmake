include(FindPackageHandleStandardArgs)

find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
  pkg_check_modules(PC_libome IMPORTED_TARGET libome QUIET)
endif()

find_path(libome_INCLUDE_DIR
  NAMES ome/ome.h
  HINTS
    ${PC_libome_INCLUDE_DIRS}
    ${CANDIA_libome_DIR}
    ${libome_ROOT}
  PATH_SUFFIXES include)

find_library(libome_LIBRARY
  NAMES ome
  HINTS
    ${PC_libome_LIBRARY_DIRS}
    ${CANDIA_libome_DIR}
    ${libome_ROOT}
  PATH_SUFFIXES lib lib64)

if (PC_libome_VERSION)
  set(libome_VERSION "${PC_libome_VERSION}")
endif()

find_package_handle_standard_args(libome
  REQUIRED_VARS libome_INCLUDE_DIR libome_LIBRARY
  VERSION_VAR libome_VERSION)

if (libome_FOUND AND NOT TARGET libome::libome)
  add_library(libome::libome UNKNOWN IMPORTED)
  set_target_properties(libome::libome PROPERTIES
	IMPORTED_LOCATION "${libome_LIBRARY}"
	INTERFACE_INCLUDE_DIRECTORIES "${libome_INCLUDE_DIR}")

  if (PC_libome_FOUND)
	set_property(TARGET libome::libome PROPERTY
	  INTERFACE_COMPILE_OPTIONS "${PC_libome_CFLAGS_OTHER}")
	set_property(TARGET libome::libome PROPERTY
	  INTERFACE_LINK_LIBRARIES "${PC_libome_LIBRARIES_OTHER}")
  endif()
endif()


