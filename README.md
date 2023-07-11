The Linear Algebra Cookbook
===========================
[![Ubuntu](https://github.com/akleb/LACookbook/actions/workflows/ubuntu.yaml/badge.svg)](https://github.com/akleb/LACookbook/actions/workflows/ubuntu.yaml)
[![MacOS](https://github.com/akleb/LACookbook/actions/workflows/macos.yaml/badge.svg)](https://github.com/akleb/LACookbook/actions/workflows/macos.yaml)
[![Memcheck](https://github.com/akleb/LACookbook/actions/workflows/memcheck.yaml/badge.svg)](https://github.com/akleb/LACookbook/actions/workflows/memcheck.yaml)

## Description
A linear algebra package for cooking up some linear solvers.
This is intended to be an upstream package for CFD solvers that need linear algebra operations or solvers.

## Installation
This package utilizes a typical CMake build structure.
To build we create a build directory, configure, make and install
```bash
cd /path/to/build/location/
cmake /path/to/LAC/
make
make install
```
Common options that may be used for the `cmake` command:
| Flag | Default | Options | Effect |
| :---: | :---: | :---: | :---: |
| `-DCMAKE_BUILD_TYPE` | `RelWithDebInfo` | `Debug\|Release\|RelWithDebInfo\|MinSizeRel`| Specifies default CMake optimization and debug compiler flags|
| `-DCMAKE_INSTALL_PREFIX` | `/usr/local/` | `/path/to/install/location` | Specifies the directory where the headers, library and CMake config files will be stored|
| `-DBUILD_TESTS=ON` | `OFF` | `ON\|OFF` | Set to `ON` to build tests. Tests suite can be run using `make test` after building|

## Usage
After installing, downstream usage is straightforward. 
With a CMake build system, the `find_package` approach is usually sufficient.
```cmake
# Find the Cookbook Library
find_package(LAC REQUIRED)

# create a target that needs LAC
add_library(${MY_PROJECT} ${MY_SOURCE})

# link against the Cookbook (Note that this updates include directories as well)
target_link_libraries(${MY_PROJCET} LAC:::LAC)
```
If the Cookbook was installed globally, this approach will work out of the box.
However, if Cookbook was installed somewhere locally, CMake needs to be told where to look for it.
This is done by setting the `-DCMAKE_PREFIX_PATH` option for the downstream configure to the installation location of the Cookbook.
Multiple upstream libraries can be specified by delimiting with `;`.
