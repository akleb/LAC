cmake_minimum_required(VERSION 3.10)

project(LAC VERSION 0.0.0)

enable_language(CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(MPI REQUIRED)

# get directory variables
include(LACDirs)

# set include directories
include_directories(${INCDIRS})

IF (NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
     "Build types are: Debug Release RelWithDebInfo MinSizeRel."
     FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)
IF (NOT BUILD_TESTS)
    SET(BUILD_TESTS CACHE BOOL
     "Build testing."
     FORCE)
ENDIF (NOT BUILD_TESTS)
IF (NOT MEMCHECK)
    SET(MEMCHECK CACHE BOOL
     "Enable nonzero error codes from memcheck."
     FORCE)
ENDIF (NOT MEMCHECK)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# --------------
#  Unit Testing
# --------------
IF (BUILD_TESTS)
    message(STATUS "Testing enabled.")
    IF(MEMCHECK)
        message(STATUS "Setting valgrind error code to 1.")
        SET(MEMORYCHECK_COMMAND_OPTIONS "--error-exitcode=1 --leak-check=full")
    ENDIF(MEMCHECK)
    include(CTest)
    include(add_unit_test)
    include(add_unit_test_mpi)
    # Copy all the input files to the bin for tests
    file(COPY ${INPUT_FILE_DIR} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    add_subdirectory(test)
ENDIF (BUILD_TESTS)

# --------------------
#  Libraries to Build
# --------------------
add_compile_options(-Wall -Werror)
add_subdirectory(src)

