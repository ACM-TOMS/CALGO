# Name the project you can choose any name you want here
PROJECT(lap) 

# Check requirements
CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
FIND_PACKAGE(CUDA)

# set the compiler flags
SET(CMAKE_CXX_COMPILER g++-7)
SET(CMAKE_C_COMPILER gcc-7)

# Release Mode
if (${CUDA_FOUND})
  SET(CUDA_PROPAGATE_HOST_FLAGS "OFF")
endif (${CUDA_FOUND})
SET(CMAKE_CXX_FLAGS "-Wall -O3 -Wfatal-errors -fstrict-aliasing -m64 -std=c++11 -fopenmp -flto -march=native -mfpmath=sse")
if (${CUDA_FOUND})
  SET(CUDA_NVCC_FLAGS "-arch=sm_52 --ptxas-options=-v -Xcompiler -Wall -O3 -std=c++11 -Xcompiler -fopenmp -Xcompiler -march=native -Xcompiler -mfpmath=sse --expt-extended-lambda")
endif (${CUDA_FOUND})

# Define where your executables should be put
SET(EXECUTABLE_OUTPUT_PATH ${MyProject_BINARY_DIR})

# Tell your source files here 
SET(test_cpu_SRCS
	../test/test_cpu.cpp
)
if (${CUDA_FOUND})
  SET(test_gpu_SRCS
          ../test/test_gpu.cu
)
endif (${CUDA_FOUND})

# Include directories
include_directories(../)
if (${CUDA_FOUND})
  include_directories(../cub-1.8.0/)
endif (${CUDA_FOUND})

ADD_EXECUTABLE(test_cpu ${test_cpu_SRCS})
# second build for counting evaluations
ADD_EXECUTABLE(test_cpu_evaluated ${test_cpu_SRCS})
TARGET_COMPILE_DEFINITIONS(test_cpu_evaluated PUBLIC LAP_DISPLAY_EVALUATED)

if (${CUDA_FOUND})
  CUDA_ADD_EXECUTABLE(test_gpu ${test_gpu_SRCS})
  # second build for counting evaluations
  CUDA_ADD_EXECUTABLE(test_gpu_evaluated ${test_gpu_SRCS})
  target_compile_definitions(test_gpu_evaluated PUBLIC LAP_DISPLAY_EVALUATED)
endif (${CUDA_FOUND})


