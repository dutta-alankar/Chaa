# Chaa — CMake build
#
#   cmake -B build                 # configure (HDF5 on by default)
#   cmake --build build            # -> build/bin/chaa
#   ctest --test-dir build         # run the validated test suite
#
# Options:
#   -DCHAA_HDF5=OFF                build without the HDF5 dependency
#   -DCHAA_NG=<n>                  ghost-cell layers (compile_params.chpl)
#   -DCHAA_CHPL_FLAGS="..."        extra chpl flags (default: --fast)
#   -DHDF5_ROOT=<prefix>           HDF5 install prefix hint

cmake_minimum_required(VERSION 3.16)
project(Chaa LANGUAGES NONE)

find_program(CHPL_EXECUTABLE chpl REQUIRED
             DOC "Chapel compiler (https://chapel-lang.org)")

option(CHAA_HDF5 "Enable HDF5 output (links -lhdf5)" ON)
set(CHAA_NG "3" CACHE STRING "ghost-cell layers (3 needed for ppm/wenoz)")
set(CHAA_NSCAL "1" CACHE STRING "passive tracer fields")
set(CHAA_CHPL_FLAGS "--fast" CACHE STRING "chpl optimization/extra flags")

file(GLOB CHAA_SOURCES
     ${CMAKE_CURRENT_SOURCE_DIR}/src/*.chpl
     ${CMAKE_CURRENT_SOURCE_DIR}/src/problems/*.chpl)
set(CHAA_H5_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/src/chaa_h5.h)

separate_arguments(CHPL_FLAG_LIST UNIX_COMMAND "${CHAA_CHPL_FLAGS}")
list(APPEND CHPL_FLAG_LIST -sNG=${CHAA_NG} -sNSCAL=${CHAA_NSCAL})

if(CHAA_HDF5)
  find_library(HDF5_LIBRARY
               NAMES hdf5
               HINTS ${HDF5_ROOT} $ENV{HDF5_ROOT} $ENV{HDF5_DIR}
               PATH_SUFFIXES lib
               PATHS /opt/homebrew/opt/hdf5/lib
                     /usr/local/opt/hdf5/lib
                     /usr/lib/x86_64-linux-gnu/hdf5/serial
                     /usr/lib/aarch64-linux-gnu/hdf5/serial)
  if(NOT HDF5_LIBRARY)
    message(FATAL_ERROR
            "libhdf5 not found; set -DHDF5_ROOT=<prefix> or -DCHAA_HDF5=OFF")
  endif()
  get_filename_component(HDF5_LIBDIR ${HDF5_LIBRARY} DIRECTORY)
  list(APPEND CHPL_FLAG_LIST -shdf5Enabled=true -lhdf5 -L${HDF5_LIBDIR})
  message(STATUS "HDF5: ${HDF5_LIBRARY}")
else()
  message(STATUS "HDF5 output disabled")
endif()

set(CHAA_BIN ${CMAKE_BINARY_DIR}/bin/chaa)

add_custom_command(
  OUTPUT ${CHAA_BIN}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
  COMMAND ${CHPL_EXECUTABLE} ${CHPL_FLAG_LIST}
          --main-module Chaa ${CHAA_SOURCES} ${CHAA_H5_HEADER}
          -o ${CHAA_BIN}
  DEPENDS ${CHAA_SOURCES} ${CHAA_H5_HEADER}
  COMMENT "Compiling Chaa with chpl"
  VERBATIM)

add_custom_target(chaa ALL DEPENDS ${CHAA_BIN})

# ------------------------- test suite (ctest) -------------------------
enable_testing()
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/tests/cases.conf CASE_LINES)
foreach(line IN LISTS CASE_LINES)
  if(line MATCHES "^[ \t]*#" OR line MATCHES "^[ \t]*$")
    continue()
  endif()
  string(REGEX REPLACE "\\|.*$" "" cname "${line}")
  string(STRIP "${cname}" cname)
  add_test(NAME ${cname}
           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_case.sh ${cname})
  set_tests_properties(${cname} PROPERTIES
    ENVIRONMENT "CHAA_BIN=${CHAA_BIN};TEST_OUT=${CMAKE_BINARY_DIR}/test-output")
endforeach()
