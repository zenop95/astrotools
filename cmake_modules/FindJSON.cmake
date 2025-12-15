# - Locate JSON library (header only)
# Defines:
#
#  JSON_FOUND
#  JSON_INCLUDE_DIR
#  JSON_INCLUDE_DIRS (not cached)
#  JSON_LIBRARY
#  JSON_LIBRARIES (not cached)

find_path(JSON_INCLUDE_DIR json/json.h
          PATH_SUFFIXES "include" "include/jsoncpp"
          DOC "JSON include directory")
mark_as_advanced(JSON_INCLUDE_DIR)

find_library(JSON_LIBRARY NAMES jsoncpp
             DOC "JSON libraries")
mark_as_advanced(JSON_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set JSON_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(JSON
                                  FOUND_VAR JSON_FOUND
                                  REQUIRED_VARS JSON_INCLUDE_DIR JSON_LIBRARY)

mark_as_advanced(JSON_FOUND)

if(JSON_FOUND)
    set(JSON_INCLUDE_DIRS ${JSON_INCLUDE_DIR})
    set(DACE_LIBRARIES ${DACE_LIBRARY})
endif()
