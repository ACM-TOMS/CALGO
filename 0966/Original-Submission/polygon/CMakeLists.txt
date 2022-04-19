# This is the CMake script for compiling a CGAL application.
project( Polygons_Basic ) 
set(CMAKE_CXX_FLAGS " -frounding-math -D_REENTRANT -DIL_STD -O3 -DNDEBUG")

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

#set(CMAKE_VERBOSE_MAKEFILE true)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS )
include( ${CGAL_USE_FILE} )

if ( CGAL_FOUND )

  # Add the libraries
  add_library(PolygonExt SHARED PolygonExt.C)

  add_library(PolygonWithHolesExt SHARED PolygonWithHolesExt.C)
  target_link_libraries(PolygonWithHolesExt PolygonExt)

  add_library(polAlgorithms SHARED polAlgorithms.C)

else()

  message(STATUS "NOTICE: This demo requires CGAL  and will not be compiled.")

endif()
