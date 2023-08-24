set(LAC_SRC     ${LAC_SRC_DIR}/lac_Default.cpp
                ${LAC_SRC_DIR}/lac_Direct.cpp
                ${LAC_SRC_DIR}/lac_Error.cpp
                ${LAC_SRC_DIR}/lac_GMRES.cpp
                ${LAC_SRC_DIR}/lac_MatrixMath.cpp
                ${LAC_SRC_DIR}/lac_MPI.cpp
                ${LAC_SRC_DIR}/lac_Norms.cpp
                ${LAC_SRC_DIR}/lac_Version.cpp)

execute_process(COMMAND git describe --always --tags --dirty OUTPUT_VARIABLE GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "The Git version is ${GIT_VERSION}")
add_definitions(-DLAC_GIT_VERSION="${GIT_VERSION}")

# ----------------
#  Library Target
# ----------------
# allows us to get global install location by default and updated by
# CMAKE_INSTALL_PREFIX option. This also allows package managers to specify
# install locations
include(GNUInstallDirs)

# build static library from source files
add_library(LAC STATIC ${LAC_SRC})

target_include_directories(LAC PRIVATE ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(LAC ${MPI_CXX_LIBRARIES})

# We want cmake to be able to pickup on the include file locations in downstream
# projects, so we want PUBLIC includes. When it is built locally the includes
# are found relative to the project directory structure. When it is installed
# and being used in a downstream package, we want to include from the install
# location
target_include_directories(LAC PUBLIC 
  $<BUILD_INTERFACE:${LAC_INC_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# We want to export our library target, tracks where everything is installed
install(TARGETS LAC EXPORT LAC_Exports
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# This is going to write the cmake version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_BINARY_DIR}/LAC/LACConfigVersion.cmake"
  VERSION ${CMAKE_PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Exports the target information with namespace LAC:: to the LACTargets.cmake
# file with is still in the local build directory, it gets installed during an
# install command. This file is included by the LACConfig, which is what CMake
# looks for when it looks for a package
export(TARGETS LAC
  FILE "${CMAKE_BINARY_DIR}/LAC/LACTargets.cmake"
  NAMESPACE LAC::
)
# Copies the LACConfig into the build directory, so that all the .cmakes can be
# installed together
configure_file(${LAC_TOPDIR}/cmake/LACConfig.cmake
  "${CMAKE_BINARY_DIR}/LAC/LACConfig.cmake"
  COPYONLY
)

# A lib/cmake/LAC is a directory structure cmake will search for, we store it in
# this variable
set(ConfigPackageLocation ${CMAKE_INSTALL_LIBDIR}/cmake/LAC)

# Install the exported targets to the install location rather than build
# directory (could be the same if the user sets CMAKE_INSTALL_PREFIX to the
# build directory.
install(EXPORT LAC_Exports
  FILE
    LACTargets.cmake
  NAMESPACE
    LAC::
  DESTINATION
    ${ConfigPackageLocation}
)
# Installs the Version and Default Config file to their install location
install(
  FILES
    ${CMAKE_BINARY_DIR}/LAC/LACConfig.cmake
    ${CMAKE_BINARY_DIR}/LAC/LACConfigVersion.cmake
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)

# ------------------
#  Install Headers
# ------------------
install(
  DIRECTORY ${LAC_INC_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Devel
)

