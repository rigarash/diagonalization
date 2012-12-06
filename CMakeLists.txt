#  Copyright Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp> 2012.
#  Distributed under the Boost Software License, Version 1.0.
#      (See http://www.boost.org/LICENSE_1_0.txt)
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(diag-parapack NONE)

######################################################################
# options
######################################################################

option(DIAG_USE_MKL_PARALLEL "Use parallelized routine of Intel MKL" ON)
if(DIAG_USE_MKL_PARALLEL)
  add_definitions(-DDIAG_USE_MKL_PARALLEL)
endif(DIAG_USE_MKL_PARALLEL)

######################################################################
# find ALPS Library
######################################################################
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

######################################################################
# applications
######################################################################
# enable C and C++ compilers
enable_language(CXX)

# diagonalization convenience library
add_library(diag fulldiag.cc sparsediag.cc)

# diag-parapack
add_executable(diag-parapack main.cc)
target_link_libraries(diag-parapack diag ${ALPS_LIBRARIES})
#add_alps_test(diag-parapack)

# diag-parapack-ns
add_executable(diag-parapack-ns diag_ns.cc)
target_link_libraries(diag-parapack-ns diag ${ALPS_LIBRARIES})
#add_alps_test(diag-parapack-ns)
