include(FindPackageHandleStandardArgs)

find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
  pkg_check_modules(LIBOME libome QUIET)
  if (LIBOME_FOUND)
	set(LIBOME_INCLUDE_DIRS ${LIBOME_INCLUDE_DIRS} CACHE STRING INTERNAL)
	set(LIBOME_LIBRARIES ${LIBOME_LIBRARIES} CACHE STRING INTERNAL)
  endif()
endif()
if (NOT PkgConfig_FOUND OR NOT LIBOME_FOUND)
  find_path(LIBOME_INCLUDE_DIR
	include/ome/ome.h
	HINTS ${CANDIA_LIBOME_DIR})
  if (LIBOME_INCLUDE_DIR)
	set(LIBOME_INCLUDE_DIRS ${LIBOME_INCLUDE_DIR}/include CACHE STRING INTERNAL)
  else()
	message("-- candia-v2: failed to find include dir ${CANDIA_LIBOME_DIR}/include/ome/ome.h")
  endif()
  
  find_library(LIBOME_LIBRARY
	NAMES ome
	HINTS ${CANDIA_LIBOME_DIR}/lib ${CANDIA_LIBOME_DIR}/lib64
	NO_DEFAULT_PATHS)
  if (LIBOME_LIBRARY)
	set(LIBOME_LIBRARIES ${LIBOME_LIBRARY} CACHE STRING INTERNAL)
  else()
	message("-- candia-v2: failed to find lib ${CANDIA_LIBOME_DIR}/lib(64)/libome.a")
  endif()
endif()

find_package_handle_standard_args(LIBOME REQUIRED_VARS LIBOME_INCLUDE_DIRS LIBOME_LIBRARIES)
