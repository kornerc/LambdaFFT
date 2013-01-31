project(lfft)

cmake_minimum_required(VERSION 2.6)

set(LFFT_VERSION_MAJOR 0)
set(LFFT_VERSION_MINOR 1)
set(LFFT_VERSION_PATCH 0)
set(LFFT_VERSION ${LFFT_VERSION_MAJOR}.${LFFT_VERSION_MINOR}.${LFFT_VERSION_PATCH})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

option(BUILD_UNIT_TESTS "Build the unit tests." OFF)

add_subdirectory(src)

if(BUILD_UNIT_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif(BUILD_UNIT_TESTS)

