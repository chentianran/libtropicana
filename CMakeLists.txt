cmake_minimum_required (VERSION 2.6)


project (trop CXX)

#set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DBOOST_NOINLINE='__attribute__ ((noinline))'")

find_package (Boost 1.44.0 REQUIRED COMPONENTS system timer)

find_package (Eigen3 REQUIRED)
#find_package (CUDA REQUIRED)

#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -DBOOST_NOINLINE='__attribute__ ((noinline))')

include_directories ("${Boost_INCLUDE_DIRS}")
include_directories ("${EIGEN3_INCLUDE_DIR}")
include_directories ("${CMAKE_CURRENT_LIST_DIR}")

add_executable (trop
    inverse.hh
    lookup.hh
    state.hh state.cc
    main.cc)

target_link_libraries (trop ${Boost_LIBRARIES})
