cmake_minimum_required(VERSION 3.10.0)
project(uncertain LANGUAGES CXX C VERSION 0.1.0)

#=============================================================================
# Augment CMake with our custom scripts
#=============================================================================
set(EXTRA_MODULES_DIR ${CMAKE_CURRENT_LIST_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${EXTRA_MODULES_DIR})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})

#=============================================================================
# ensure build type first
#=============================================================================

include(EnsureBuildType)

#=============================================================================
# Conan
#=============================================================================
# gersemi: off
if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
  message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
  file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/0.18.1/conan.cmake"
    "${CMAKE_BINARY_DIR}/conan.cmake"
    TLS_VERIFY ON)
endif()

set(CONAN_PROFILE "default" CACHE STRING "Name of conan profile to use, uses default by default")
set(CONAN "AUTO" CACHE STRING "conan options AUTO (conan must be in path), MANUAL (expects conanbuildinfo.cmake in build directory) or DISABLE")
if(${CONAN} MATCHES "AUTO")
  include(${CMAKE_BINARY_DIR}/conan.cmake)
  conan_cmake_run(CONANFILE conanfile.txt
    PROFILE ${CONAN_PROFILE}
    BASIC_SETUP NO_OUTPUT_DIRS KEEP_RPATHS
    BUILD_TYPE "None"
    BUILD outdated)
elseif(${CONAN} MATCHES "MANUAL")
  if(EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
    include("${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
    conan_basic_setup(NO_OUTPUT_DIRS SKIP_RPATH)
  else()
    message(FATAL_ERROR "CONAN set to MANUAL but no file named conanbuildinfo.cmake found in build directory")
  endif()
elseif(NOT ${CONAN} MATCHES "DISABLE")
  message(FATAL_ERROR "Unrecognised option for CONAN (${CONAN}), use AUTO, MANUAL or DISABLE")
endif()
# gersemi: on
#=============================================================================
# General configuration
#=============================================================================

include(CoverageReports)
include(OutputDirConfig)
include(PackageConfig)
include(BuildInfo)
#include(DoxygenUtilities)

#=============================================================================
# compiler and library configuration
#=============================================================================
include(CompilerConfig)

add_subdirectory(source)
add_subdirectory(examples)

#=============================================================================
# install the examples directory to the documentation directory
#=============================================================================
#install(DIRECTORY examples
#        DESTINATION ${CMAKE_INSTALL_DOCDIR})

#=============================================================================
# unit tests if GTest if present
#=============================================================================
find_package(Threads REQUIRED)
find_package(GTest REQUIRED)
enable_testing()
add_subdirectory(tests)

add_custom_target(everything DEPENDS uncertain unit_tests examples)
