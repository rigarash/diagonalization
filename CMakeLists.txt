cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(diag-parapack NONE)

# find ALPS Library
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

# enable C and C++ compilers
enable_language(CXX)

# diag-parapack
add_executable(diag-parapack main.cc fulldiag.cc sparsediag.cc)
target_link_libraries(diag-parapack ${ALPS_LIBRARIES})
add_alps_test(diag-parapack)

# diag-parapack-ns
add_executable(diag-parapack-ns diag_ns.cc fulldiag.cc sparsediag.cc)
target_link_libraries(diag-parapack-ns ${ALPS_LIBRARIES})
add_alps_test(diag-parapack-ns)
