# - Locate DACE library
# Defines:
#
#  DACE_FOUND
#  DACE_INCLUDE_DIR
#  DACE_INCLUDE_DIRS (not cached)
#  DACE_LIBRARY
#  DACE_LIBRARIES (not cached)

find_path(DACE_INCLUDE_DIR /dace/dace.h
          HINTS $ENV{DACE_ROOT_DIR}/include ${DACE_ROOT_DIR}/include
          DOC "DACE include directory")
mark_as_advanced(DACE_INCLUDE_DIR)

find_library(DACE_LIBRARY NAMES dace
             HINTS $ENV{DACE_ROOT_DIR}/lib ${DACE_ROOT_DIR}/lib
             DOC "DACE libraries")
mark_as_advanced(DACE_LIBRARY DACE_CBLAS_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set DACE_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DACE
                                  FOUND_VAR DACE_FOUND
                                  REQUIRED_VARS DACE_INCLUDE_DIR DACE_LIBRARY)

mark_as_advanced(DACE_FOUND)

if(DACE_FOUND)
    set(DACE_INCLUDE_DIRS ${DACE_INCLUDE_DIR})
    set(DACE_LIBRARIES ${DACE_LIBRARY})
endif()
