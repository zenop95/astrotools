#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "cfg::cfg" for configuration "RelWithDebInfo"
set_property(TARGET cfg::cfg APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(cfg::cfg PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libcfg.so.0.1.0"
  IMPORTED_SONAME_RELWITHDEBINFO "libcfg.so.0.1.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS cfg::cfg )
list(APPEND _IMPORT_CHECK_FILES_FOR_cfg::cfg "${_IMPORT_PREFIX}/lib/libcfg.so.0.1.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
