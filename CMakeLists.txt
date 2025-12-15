cmake_minimum_required (VERSION 3.10.0 FATAL_ERROR)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)

# Define project
project(CDTI
        VERSION 0.1.0
        DESCRIPTION "CDTI project software"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)
include(cdtiUtils)

# Compiler options
add_compile_options("-Wall")
add_compile_options("-Wextra")
add_compile_options("-Wpedantic")
add_compile_options("-Wfloat-equal")
add_compile_options("-Wnull-dereference")
add_compile_options("-fPIC")
set(CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

# Set the type of build
set(default_build_type "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for each cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
               "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

# Project organized by folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

message(STATUS "CMAKE_SYSTEM_NAME          = ${CMAKE_SYSTEM_NAME}"         )
message(STATUS "CMAKE_CXX_COMPILER_ID      = ${CMAKE_CXX_COMPILER_ID}"     )
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_COMPILER_IS_GNUCXX   = ${CMAKE_COMPILER_IS_GNUCXX}"  )
message(STATUS "BUILD TYPE                 = ${CMAKE_BUILD_TYPE}"          )

# Spice libraries
find_package(CSPICE)
message(STATUS "CSPICE_FOUND               = ${CSPICE_FOUND}" )
include_directories(${CSPICE_INCLUDE_DIR})
add_definitions(-DUSE_SPICE)

# General testing information
find_package(CPPUNIT 1.14.0)
message(STATUS "CPPUNIT_FOUND              = ${CPPUNIT_FOUND}" )

# DLib libraries
find_package(dlib REQUIRED)
message(STATUS "DLIB_FOUND                 = ${dlib_FOUND}" )
message(STATUS "DLIB_VERSION               = ${dlib_VERSION}" )

# Eigen libraries
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
message(STATUS "EIGEN_FOUND                = ${Eigen3_FOUND}" )
message(STATUS "EIGEN_VERSION              = ${Eigen3_VERSION}" )

find_package(BLAS REQUIRED)
message(STATUS "BLAS_FOUND                 = ${BLAS_FOUND}" )

# JSON for file parsing
find_package(JSON REQUIRED)
message(STATUS "JSON_FOUND                 = ${JSON_FOUND}" )

# Mathematical libraries
find_package(DACE 2.0.0 REQUIRED)
message(STATUS "DACE_FOUND                 = ${DACE_FOUND}" )



# Add subdirectories:
add_subdirectory(modules)
add_subdirectory(tools)

# Add test subdirectory and enable tests
include(CTest)
if(BUILD_TESTING)
    if(${CPPUNIT_FOUND})
        add_subdirectory(tests)
    else()
        message(SEND_ERROR "Cannot build unit tests, required dependency CPPunit not found.")
    endif()
endif()
