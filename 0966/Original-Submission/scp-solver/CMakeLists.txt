# This is the CMake script for compiling a CGAL application.
project( Xpress_PLI_Solver ) 
set(CMAKE_CXX_FLAGS " -frounding-math -D_REENTRANT -DIL_STD -O3 -DNDEBUG")

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

#set(CMAKE_VERBOSE_MAKEFILE true)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL)
include( ${CGAL_USE_FILE} )

set(XPRESS_PATH "${XPRESS_PATH}")
set(GLPK_PATH "${GLPK_PATH}")

if(XPRESS_PATH STREQUAL "")
    MESSAGE("XPRESS WILL NOT BE USED!")
else()
    MESSAGE("XPRESS_PATH is ${XPRESS_PATH}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXPRESS")
    
    include_directories(${XPRESS_PATH}/include)
    link_directories(${XPRESS_PATH}/lib)

    add_library(SolverPLIXpress SHARED SolverPLIXpress.C)
    target_link_libraries(SolverPLIXpress xprs)
endif()

if(GLPK_PATH STREQUAL "")
    MESSAGE("GLPK WILL NOT BE USED!")
else()
    MESSAGE("GLPK_PATH is ${XPRESS_PATH}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGLPK")
    
    include_directories(${GLPK_PATH}/include)
    link_directories(${GLPK_PATH}/lib)

    add_library(SolverPLIGlpk SHARED SolverPLIGlpk.C)
    target_link_libraries(SolverPLIGlpk glpk)
endif()

# Add the libraries
add_library(Lagrangian SHARED Lagrangian.C)

