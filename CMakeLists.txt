cmake_minimum_required(VERSION 3.6)

project(transflow)

set(TRANSFLOW_VERSION_MAJOR 0)
set(TRANSFLOW_VERSION_MINOR 1)
set(TRANSFLOW_VERSION_PATCH 0)
set(TRANSFLOW_VERSION "${TRANSFLOW_VERSION_MAJOR}.${TRANSFLOW_VERSION_MINOR}.${TRANSFLOW_VERSION_PATCH}")

set(PROJ_NAME "\"TransFlow\"") # PROJECT_NAME is a predefined variable
set(PROJECT_DESCRIPTION "\"A transient pipeline flow simulation library\"")

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

configure_file(src/config.hpp.in
    ${CMAKE_BINARY_DIR}/src/config.hpp
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# don't use extended variants of language dialects (like -std=gnu++11 instead of -std=c++11)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(BUILD_TESTS "Build test executable" ON)
option(BUILD_DOCS "Generate documentation" ON)
option(BUILD_EXAMPLES "Build examples" ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as no build type was specified")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type (Debug/Release)" FORCE)
endif (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /WX")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
endif(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

set(LIBRARY_SOURCES
    src/pipeline.cpp
    src/equationofstate/equationofstatebase.cpp
    src/equationofstate/bwrs.cpp
    src/equationofstate/gerg04.cpp
    src/equationofstate/equationofstate.cpp
    src/heattransfer/utils.cpp
    src/heattransfer/material.cpp
    src/heattransfer/burialmedium.cpp
    src/heattransfer/ambientfluid.cpp
    src/heattransfer/pipewall.cpp
    src/heattransfer/radial.cpp
    src/heattransfer/heattransfer.cpp
    src/heattransfer/heattransferbase.cpp
    src/heattransfer/heattransferstate.cpp
    src/heattransfer/unsteady.cpp
    src/heattransfer/steadystate.cpp
    src/utilities/utilities.cpp
    src/utilities/physics.cpp
    src/utilities/numerics.cpp
    src/utilities/linearinterpolator.cpp
    src/solver/discretizer/discretizer.cpp
    src/solver/discretizer/enthalpy.cpp
    src/solver/discretizer/internalenergy.cpp
    src/solver/matrixequation.cpp
    src/solver/governingequationsolver.cpp
    src/advection/batchtracking.cpp
    src/advection/batchtrackingstate.cpp
    src/solver/solver.cpp
    src/physics.cpp
    src/simulator.cpp
    src/solver/boundaryconditions.cpp
    src/composition.cpp
    src/timeseries.cpp
    src/sampler.cpp
)

add_library(${PROJECT_NAME} ${LIBRARY_SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/src)

find_package(Armadillo REQUIRED)
target_include_directories(${PROJECT_NAME} PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${ARMADILLO_LIBRARIES})
message("ARMADILLO_VERSION_MAJOR: " ${ARMADILLO_VERSION_MAJOR})

find_package(HDF5 REQUIRED)
target_include_directories(${PROJECT_NAME} PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES})
message("HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
message("HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")

# had to adjust armadillo_bits/config.hpp, and comment out the line
#     #define ARMA_SUPERLU_INCLUDE_DIR /
# to get SuperLU to work properly
find_package(SuperLU REQUIRED)
#target_include_directories(${PROJECT_NAME} PRIVATE ${SUPERLU_INCLUDE_DIRS})  # not used
target_link_libraries(${PROJECT_NAME} ${SUPERLU_LIBRARIES})
message("SUPERLU_LIBRARIES: ${SUPERLU_LIBRARIES}")
message("SUPERLU_INCLUDE_DIRS: ${SUPERLU_INCLUDE_DIRS}")

# set up arma superlu include path
add_compile_definitions(ARMA_SUPERLU_INCLUDE_DIR=${SUPERLU_INCLUDE_DIRS}/)  # note trailing slash!
message("ARMA_SUPERLU_INCLUDE_DIR (incl. trailing slash): ${SUPERLU_INCLUDE_DIRS}/")

find_package(LAPACK REQUIRED)
target_link_libraries(${PROJECT_NAME} ${LAPACK_LIBRARIES})

target_link_libraries(${PROJECT_NAME} stdc++fs)

target_compile_definitions(${PROJECT_NAME} PUBLIC
    ARMA_DONT_USE_WRAPPER
    ARMA_USE_SUPERLU
    ARMA_USE_HDF5
    ARMA_USE_LAPACK
    ARMA_USE_BLAS
    ARMA_USE_ARPACK
    ARMA_PRINT_ERRORS
)

# set up resource folder
FILE(COPY ${CMAKE_SOURCE_DIR}/res DESTINATION ${CMAKE_BINARY_DIR})
set(TRANSFLOW_RESOURCE_PATH \"${CMAKE_BINARY_DIR}/res/\")
target_compile_definitions(${PROJECT_NAME} PUBLIC TRANSFLOW_RESOURCE_PATH=${TRANSFLOW_RESOURCE_PATH})
message("TRANSFLOW_RESOURCE_PATH: ${TRANSFLOW_RESOURCE_PATH}")

if (BUILD_TESTS)
    include(CTest)
    enable_testing()
    add_subdirectory(test)
endif (BUILD_TESTS)

if (BUILD_DOCS)
    add_subdirectory(doc)
endif (BUILD_DOCS)

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif (BUILD_EXAMPLES)
