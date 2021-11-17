cmake_minimum_required(VERSION 3.0)
project(FEA_Mesher2D)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")

#Set the location for your runtime directory, this is where the executable will be built to
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")

set(SOURCE_FILES Application.cpp Edge.cpp Ray.cpp GeoPrimitives.cpp BoundaryLayerMesh.cpp InviscidRegionMesh.cpp MeshGenerator.cpp Vertex.cpp BoundaryLayerSubdomain.cpp InviscidRegionSubdomain.cpp MeshingManager.cpp ADT2D.cpp
ADT2DElement.cpp ADT2DExtent.cpp)

set(SOURCE_FILES ${SOURCE_FILES} main.cpp)
add_library(objects OBJECT ${SOURCE_FILES})
add_executable(FEA_Mesher2D $<TARGET_OBJECTS:objects> triangle/triangle.o)
target_link_libraries(FEA_Mesher2D ${MPI_LIBRARIES})
