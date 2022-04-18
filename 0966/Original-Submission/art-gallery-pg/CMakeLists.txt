# This is the CMake script for compiling a CGAL application.
project( Art_Gallery_Solver ) 
set(CMAKE_CXX_FLAGS " -frounding-math -D_REENTRANT -DIL_STD -O3 -DNDEBUG")

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

#set(CMAKE_VERBOSE_MAKEFILE true)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

set(XPRESS_PATH "${XPRESS_PATH}")
set(GLPK_PATH "${GLPK_PATH}")

if(XPRESS_PATH STREQUAL "")
    MESSAGE("XPRESS WILL NOT BE USED!")
else()
    MESSAGE("XPRESS_PATH is ${XPRESS_PATH}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXPRESS")
    
    include_directories(${XPRESS_PATH}/include)
    link_directories(${XPRESS_PATH}/lib)
endif()

if(GLPK_PATH STREQUAL "")
    MESSAGE("GLPK WILL NOT BE USED!")
else()
    MESSAGE("GLPK_PATH is ${GLPK_PATH}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGLPK")
    
    include_directories(${GLPK_PATH}/include)
    link_directories(${GLPK_PATH}/lib)
endif()

find_package(CGAL COMPONENTS)
include( ${CGAL_USE_FILE} )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../polygon)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../scp-solver)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../grid)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../pre-solver)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../polygon)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../scp-solver)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../grid)
link_directories(${CMAKE_CURRENT_SOURCE_DIR})
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../pre-solver)
  
add_executable(artGallerySolver artGallerySolver.C AuxGallery.C)

target_link_libraries( artGallerySolver
    m Lagrangian PreSolver PolygonExt PolygonWithHolesExt polAlgorithms 
    IGrid AVPLightGrid
    )

