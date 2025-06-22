# - Try to find Eigen
# Defines:
#
#  EIGEN_FOUND
#  EIGEN_INCLUDE_DIR
#  EIGEN_INCLUDE_DIRS (not cached)

find_path(EIGEN_INCLUDE_DIR Eigen/Core
          HINTS $ENV{EIGEN_ROOT_DIR}/include/eigen3 ${EIGEN_ROOT_DIR}/include/eigen3)

mark_as_advanced(EIGEN_INCLUDE_DIR)

# handle the QUIETLY and REQUIRED arguments and set LHAPDF_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Eigen DEFAULT_MSG EIGEN_INCLUDE_DIR)

set(EIGEN_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR})

mark_as_advanced(EIGEN_FOUND)
