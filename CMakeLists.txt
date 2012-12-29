#  Copyright Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp> 2012.
#  Distributed under the Boost Software License, Version 1.0.
#      (See http://www.boost.org/LICENSE_1_0.txt)

# The minimum version is stricter than ALPS itself,
# in order to load the correct FindBoost.cmake when testing.
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(diag-parapack NONE)

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

# simple
add_executable(simple simple.cc simple_matrix_worker.cc)
target_link_libraries(simple ${ALPS_LIBRARIES})
add_alps_test(simple)

######################################################################
# unit test (using Boost.Test)
######################################################################

# Check and build Boost.Test Unit Test Framework
# since ALPS does not build it by default.
if (NOT ALPS_HAS_BOOST_TEST)
  find_package(Boost 1.48.0 COMPONENTS unit_test_framework)
endif (NOT ALPS_HAS_BOOST_TEST)

if (ALPS_HAS_BOOST_TEST OR Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  add_definitions(-DBOOST_TEST_DYN_LINK)
  add_executable(test_types test_types.cc)
  target_link_libraries(test_types ${ALPS_LIBRARIES} ${Boost_LIBRARIES})
  enable_testing()
  add_test(test_types test_types)
endif (ALPS_HAS_BOOST_TEST OR Boost_FOUND)
