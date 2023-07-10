set(LAC_SRC     ${LAC_DIR}/lac_Default.cpp
                ${LAC_DIR}/lac_Error.cpp
                ${LAC_DIR}/lac_LinearAlgebra.cpp
                ${LAC_DIR}/lac_MPI.cpp
                ${LAC_DIR}/lac_Version.cpp)

execute_process(COMMAND git describe --always --tags --dirty OUTPUT_VARIABLE GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "The Git version is ${GIT_VERSION}")
add_definitions(-DLAC_GIT_VERSION="${GIT_VERSION}")

# ----------------
#  Library Target
# ----------------
add_library(LACookbook ${LAC_SRC})
target_include_directories(LACookbook PRIVATE ${LAC_INC})
target_link_libraries(LACookbook ${MPI_CXX_LIBRARIES})

# ------------------
#  Install Headers
# ------------------
install(DIRECTORY ${TOPDIR}/include/ DESTINATION ${CMAKE_BINARY_DIR}/include)

