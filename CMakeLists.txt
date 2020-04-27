cmake_minimum_required (VERSION 2.8)

project (trop CXX)
enable_testing()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set (CMAKE_CXX_STANDARD 11)

option(USE_SPARSE "Use sparse matrix in the core algorithm" OFF)

if(USE_SPARSE)
    add_definitions(-DUSE_SPARSE)
endif()

if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER_EQUAL 3.12)
    find_package (Python3 COMPONENTS Development)
    include_directories ("${Python3_INCLUDE_DIRS}")
else (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER_EQUAL 3.12)
    find_package(PythonLibs 3.6 REQUIRED)
    include_directories ("${PYTHON_INCLUDE_DIRS}")
endif(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER_EQUAL 3.12)

find_package (Boost 1.63.0 REQUIRED COMPONENTS system timer unit_test_framework python3 numpy3)
find_package (Eigen3 REQUIRED)

include_directories ("${Boost_INCLUDE_DIRS}")
include_directories ("${EIGEN3_INCLUDE_DIR}")
include_directories ("${CMAKE_CURRENT_LIST_DIR}")

add_library(tropicana_dn STATIC inverse.hh lookup.hh debug.hh debug.cc state.hh state.cc cell.hh cell.cc trop.hh trop.cc)
add_library(tropicana_sp STATIC inverse.hh lookup.hh debug.hh debug.cc state.hh state.cc cell.hh cell.cc trop.hh trop.cc)
set_target_properties (tropicana_dn PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties (tropicana_sp PROPERTIES COMPILE_FLAGS -DUSE_SPARSE )

add_library(pytropicana SHARED python.cc)
set_target_properties(pytropicana PROPERTIES PREFIX "")
target_link_libraries (pytropicana tropicana_dn Boost::python3 Boost::numpy3)

target_link_libraries (tropicana_dn)
target_link_libraries (tropicana_sp)

add_executable(tropicon-dn main.cc)
target_link_libraries (tropicon-dn tropicana_dn)

add_executable(tropicon-sp main.cc)
target_link_libraries (tropicon-sp tropicana_sp)

#add_executable(pivot-perf tests/pivot_perf.cc)
#target_link_libraries (pivot-perf tropicana)

add_executable(tests-dn tests/square.cc tests/cube.cc tests/hypercube.cc tests/main.cc)
target_link_libraries (tests-dn tropicana_dn Boost::unit_test_framework)

add_executable(tests-sp tests/square.cc tests/cube.cc tests/hypercube.cc tests/main.cc)
target_link_libraries (tests-sp tropicana_sp Boost::unit_test_framework)

install (TARGETS tropicon-dn DESTINATION bin)
install (TARGETS tropicon-sp DESTINATION bin)
install (TARGETS tropicana_dn DESTINATION lib)
install (TARGETS tropicana_sp DESTINATION lib)
install (TARGETS pytropicana  DESTINATION lib)

