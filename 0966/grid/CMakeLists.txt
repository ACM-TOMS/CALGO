# This is the CMake script for compiling a CGAL application.
project( Discretization_Strategies ) 
set(CMAKE_CXX_FLAGS " -frounding-math -D_REENTRANT -DIL_STD -O3 -DNDEBUG")

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

#set(CMAKE_VERBOSE_MAKEFILE true)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS )
include( ${CGAL_USE_FILE} )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../polygon)

if ( CGAL_FOUND )

  add_library(IGrid SHARED IGrid.C)

  add_library(AVPLightGrid SHARED AVPLightGrid.C)

else()

  message(STATUS "NOTICE: This demo requires CGAL  and will not be compiled.")

endif()
