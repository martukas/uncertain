cmake_minimum_required(VERSION 3.0.0)
project(uncertain
    LANGUAGES CXX C
    VERSION 0.1.0
    )

#=============================================================================
# Augment CMake with our custom scripts
#=============================================================================
set(EXTRA_MODULES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${EXTRA_MODULES_DIR})
include(${EXTRA_MODULES_DIR}/EnsureBuildType.cmake)

#=============================================================================
# Conan
#=============================================================================
SET(CONAN_PROFILE "default" CACHE STRING "Name of conan profile to use, uses default by default")
SET(CONAN_MODE "AUTO" CACHE STRING "conan options AUTO (conan must be in path), MANUAL (expects conanbuildinfo.cmake in build directory) or DISABLE")
if(${CONAN_MODE} MATCHES "AUTO")
  include(${EXTRA_MODULES_DIR}/ConanSetup.cmake)
  conan_cmake_run(CONANFILE conanfile.txt
      PROFILE ${CONAN_PROFILE}
      BASIC_SETUP NO_OUTPUT_DIRS KEEP_RPATHS
      BUILD_TYPE "None"
      BUILD outdated)
elseif(${CONAN_MODE} MATCHES "MANUAL")
  if(EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
    include("${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
    conan_basic_setup(NO_OUTPUT_DIRS SKIP_RPATH)
  else()
    MESSAGE(FATAL_ERROR "CONAN_MODE set to MANUAL but no file named conanbuildinfo.cmake found in build directory")
  endif()
elseif(NOT ${CONAN_MODE} MATCHES "DISABLE")
  MESSAGE(FATAL_ERROR "Unrecognised option for CONAN_MODE (${CONAN_MODE}), use AUTO, MANUAL or DISABLE")
endif()

#=============================================================================
# General configuration
#=============================================================================

include(${EXTRA_MODULES_DIR}/CoverageReports.cmake)
include(${EXTRA_MODULES_DIR}/OutputDirConfig.cmake)
include(${EXTRA_MODULES_DIR}/PackageConfig.cmake)
include(${EXTRA_MODULES_DIR}/BuildInfo.cmake)
#include(DoxygenUtilities)

#=============================================================================
# compiler and library configuration
#=============================================================================
include(${EXTRA_MODULES_DIR}/CompilerConfig.cmake)

include(${EXTRA_MODULES_DIR}/ClangTidy.cmake)

#include(${EXTRA_MODULES_DIR}/Checks.cmake)

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
find_package(GTest REQUIRED)
add_subdirectory(tests)
