#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "geco::geco" for configuration "RelWithDebInfo"
set_property(TARGET geco::geco APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(geco::geco PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libgeco.so.0.1.0"
  IMPORTED_SONAME_RELWITHDEBINFO "libgeco.so.0.1.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS geco::geco )
list(APPEND _IMPORT_CHECK_FILES_FOR_geco::geco "${_IMPORT_PREFIX}/lib/libgeco.so.0.1.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
