cmake_minimum_required (VERSION 2.8)

project (trop CXX)
enable_testing()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set (CMAKE_CXX_STANDARD 11)

option(USE_SPARSE "Use sparse matrix in the core algorithm" OFF)

if(USE_SPARSE)
    add_definitions(-DUSE_SPARSE)
endif()

find_package (Boost 1.44.0 REQUIRED COMPONENTS system timer unit_test_framework)

find_package (Eigen3 REQUIRED)

include_directories ("${Boost_INCLUDE_DIRS}")
include_directories ("${EIGEN3_INCLUDE_DIR}")
include_directories ("${CMAKE_CURRENT_LIST_DIR}")

add_library(tropicana STATIC
    inverse.hh
    lookup.hh
    debug.hh debug.cc
    state.hh state.cc
    cell.hh cell.cc
    trop.hh trop.cc)

target_link_libraries (tropicana ${Boost_LIBRARIES})

add_executable(tropicon main.cc)
target_link_libraries (tropicon tropicana)

add_executable(pivot-perf tests/pivot_perf.cc)
target_link_libraries (pivot-perf tropicana)

add_executable(tests
    tests/square.cc
    tests/cube.cc
    tests/hypercube.cc
    tests/main.cc)
target_link_libraries (tests tropicana)
#add_test( square)
