#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "astro::astro" for configuration "RelWithDebInfo"
set_property(TARGET astro::astro APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(astro::astro PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libastro.so.0.1.0"
  IMPORTED_SONAME_RELWITHDEBINFO "libastro.so.0.1.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS astro::astro )
list(APPEND _IMPORT_CHECK_FILES_FOR_astro::astro "${_IMPORT_PREFIX}/lib/libastro.so.0.1.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
