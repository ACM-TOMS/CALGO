cmake_minimum_required(VERSION 3.17.5 FATAL_ERROR)
project(CP-CALS LANGUAGES CXX)
include(CheckCXXCompilerFlag)
include(cmake/utils.cmake)

option(WITH_MKL "Create MKL Library version" Off)
option(WITH_BLIS "Create BLIS Library version" Off)
option(WITH_OPENBLAS "Create OpenBLAS Library version" Off)
option(WITH_CUBLAS "Create cuBLAS Library version" Off)
option(WITH_MATLAB "Enable MATLAB support" Off)
option(WITH_TESTS "Enable Tests" Off)
option(WITH_TRAVIS "Travis mode" Off)  # There is a bug with cmake on travis, disable specifying BLAS version for now
option(WITH_DIAGNOSTICS "Set level of diagnostics during compilation (0, 1, 2) " 2)
option(WITH_EXPERIMENTS "Enable experiments that make use of external software for comparisons (CTF, PLANC)" Off)

if (WITH_DIAGNOSTICS GREATER_EQUAL 1)
    message(STATUS "Enabling verbose Makefile")
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif ()

#########################################################################################
# Find OpenMP
#########################################################################################
FIND_PACKAGE(OpenMP REQUIRED COMPONENTS CXX)
message(STATUS "${OpenMP_CXX_LIBRARIES}")

#########################################################################################
# Find CUBLAS
#########################################################################################
if (WITH_CUBLAS)
    ENABLE_LANGUAGE(CUDA)
    FIND_PACKAGE(CUDAToolkit REQUIRED)
endif ()

#########################################################################################
# Find MKL
#########################################################################################
if (WITH_MKL)
    if (NOT WITH_TRAVIS)
        set(BLA_VENDOR Intel10_64lp)
    endif ()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    # If compiler is GNU, make sure that CALS links to the GOMP version of MKL
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

        # TODO: Supposed to be a better way of doing this using list(TRANSFORM ...), but I couldn't figure it out
        set(_BLAS_LIBRARIES)
        foreach (INC ${BLAS_LIBRARIES})
            string(REGEX REPLACE ".*/libiomp5.so" "-lgomp" INC ${INC})
            string(REPLACE "intel_thread" "gnu_thread" INC ${INC})
            LIST(APPEND _BLAS_LIBRARIES ${INC})
        endforeach ()
        set(BLAS_LIBRARIES ${_BLAS_LIBRARIES})
        set(_LAPACK_LIBRARIES)
        foreach (INC ${LAPACK_LIBRARIES})
            string(REGEX REPLACE ".*/libiomp5.so" "-lgomp" INC ${INC})
            string(REPLACE "intel_thread" "gnu_thread" INC ${INC})
            LIST(APPEND _LAPACK_LIBRARIES ${INC})
        endforeach ()
        set(LAPACK_LIBRARIES ${_LAPACK_LIBRARIES})
        message(STATUS "${BLAS_LIBRARIES}")
        message(STATUS "${LAPACK_LIBRARIES}")
    endif ()

    set(BLA_VEN MKL)
endif ()

#########################################################################################
# Find BLIS
#########################################################################################
if (WITH_BLIS)
    if (NOT WITH_TRAVIS)
        set(BLA_VENDOR FLAME)
    endif ()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    set(BLA_VEN FLAME)
endif ()

#########################################################################################
# Find OpenBLAS
#########################################################################################
if (WITH_OPENBLAS)
    if (NOT WITH_TRAVIS)
        set(BLA_VENDOR OpenBLAS)
    endif ()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    set(BLA_VEN OPENBLAS)
endif ()

#########################################################################################
# Find MATLAB
#########################################################################################
if (WITH_MATLAB)
    option(MATLAB_PATH "Path to MATLAB installation")
    if (MATLAB_PATH)
        SET(Matlab_ROOT_DIR "${MATLAB_PATH}")
    else ()
        SET(Matlab_ROOT_DIR "$ENV{HOME}/opt/Matlab2018b")
    endif ()

    message("MATLAB_ROOT_DIR: " ${MATLAB_PATH})
    find_package(Matlab COMPONENTS MX_LIBRARY REQUIRED)

    set(CALS_LINK_LIBRARIES ${CALS_LINK_LIBRARIES} ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY})

    set(HAVE_MATLAB ON)

    set(BLA_VEN MATLAB)
endif ()

#########################################################################################
# Configuration
#########################################################################################
set(CMAKE_CXX_STANDARD 17)

# CP_CALS Library
set(CP_CALS CALS_${BLA_VEN})
add_library(${CP_CALS} SHARED)

target_sources(${CP_CALS}
        PRIVATE
        src/als.cpp
        src/cals.cpp
        src/tensor.cpp
        src/matrix.cpp
        src/ktensor.cpp
        src/cals_blas.cpp
        src/cuda_utils.cpp
        src/utils/utils.cpp
        src/utils/error.cpp
        src/utils/mttkrp.cpp
        src/utils/update.cpp
        src/multi_ktensor.cpp
        src/utils/line_search.cpp

        PUBLIC
        include/als.h
        include/cals.h
        include/tensor.h
        include/matrix.h
        include/ktensor.h
        include/cals_blas.h
        include/cuda_utils.h
        include/utils/utils.h
        include/utils/error.h
        include/utils/mttkrp.h
        include/utils/update.h
        include/multi_ktensor.h
        include/utils/line_search.h

        include/timer.h
        include/definitions.h
        )

target_compile_definitions(${CP_CALS} PUBLIC ${CP_CALS}=1)  # Enable the appropriate BLAS/LAPACK headers in the code

set_target_properties(${CP_CALS} PROPERTIES POSITION_INDEPENDENT_CODE ON)  # -fPIC

target_include_directories(${CP_CALS} PUBLIC include include/utils)

target_link_libraries(${CP_CALS} PUBLIC OpenMP::OpenMP_CXX)

target_compile_definitions(${CP_CALS} PUBLIC SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# Link to BLAS and LAPACK and enable detailed timing, if not in MATLAB.
if (NOT HAVE_MATLAB)
    target_compile_definitions(${CP_CALS} PUBLIC WITH_TIME=1)
    target_link_libraries(${CP_CALS} PUBLIC ${BLAS_LIBRARIES})
    target_link_libraries(${CP_CALS} PUBLIC ${LAPACK_LIBRARIES})
else ()
    target_include_directories(${CP_CALS} PRIVATE ${Matlab_INCLUDE_DIRS})
    target_compile_definitions(${CP_CALS} PUBLIC WITH_TIME=0)
    target_link_libraries(${CP_CALS} PUBLIC ${Matlab_MX_LIBRARY})
    target_link_libraries(${CP_CALS} PUBLIC ${Matlab_MEX_LIBRARY})
endif ()

################################
# Optimization Flags
################################
# TODO: Use -w and -fcompare-debug-second to eliminate Cyclops warnings during compilation
set(comp_flags -Wall)
if (CMAKE_BUILD_TYPE MATCHES Release)
    list(APPEND comp_flags -ffast-math -march=native)

    if (${BLA_VEN} STREQUAL "FLAME")
        list(APPEND comp_flags -Wno-unused-function)
    endif ()

    if (WITH_DIAGNOSTICS GREATER_EQUAL 1)
        list(APPEND comp_flags -g)
    endif ()

    if (WITH_DIAGNOSTICS GREATER_EQUAL 2)
        if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
            list(APPEND comp_flags -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize)
        endif ()
    endif ()
endif ()

add_compilation_flags(${CP_CALS} "${comp_flags}")

################################
# CUBLAS
################################

if (WITH_CUBLAS)
    set_property(TARGET ${CP_CALS} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET ${CP_CALS} PROPERTY CUDA_ARCHITECTURES 50 52 60 61 70 75 80)
    target_compile_definitions(${CP_CALS} PUBLIC CUDA_ENABLED=1)
    target_link_libraries(${CP_CALS} PUBLIC CUDA::cudart CUDA::cublas)

    set(CP_CALS_CUDA CALS_${BLA_VEN}_CUDA)
    add_library(${CP_CALS_CUDA} SHARED src/utils/khatri_rao.cu)
    set_property(TARGET ${CP_CALS_CUDA} PROPERTY POSITION_INDEPENDENT_CODE ON)  # -fPIC
    set_property(TARGET ${CP_CALS_CUDA} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET ${CP_CALS_CUDA} PROPERTY CUDA_ARCHITECTURES 50 52 60 61 70 75 80)
    target_link_libraries(${CP_CALS} PUBLIC ${CP_CALS_CUDA})
endif ()

################################
# GTest
################################
if (WITH_TESTS)
    find_package(GTest REQUIRED)
    enable_testing()

    include_directories(${GTEST_INCLUDE_DIRS})
    set(GTEST_LINK_FLAGS ${GTEST_BOTH_LIBRARIES} -lpthread)

    add_subdirectory(tests)
endif ()

################################
# Matlab
################################
IF (HAVE_MATLAB)
    add_library(cals_matlab_interface_library SHARED
            ${CMAKE_SOURCE_DIR}/matlab/matlab.h
            ${CMAKE_SOURCE_DIR}/matlab/matlab.cpp
            ${CMAKE_SOURCE_DIR}/matlab/matlab_parsing.h
            ${CMAKE_SOURCE_DIR}/matlab/matlab_parsing.cpp)
    target_include_directories(cals_matlab_interface_library PRIVATE ${Matlab_INCLUDE_DIRS})
    target_link_libraries(cals_matlab_interface_library CALS_MATLAB ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY})
    set(CALS_MATLAB_LINK_LIBS cals_matlab_interface_library CALS_MATLAB)

    # Newer versions of CMake have an API version argument to MATLAB_ADD_MEX
    # which needs to be set to the right value since it defaults to an incorrect
    # value
    set(CALS_MATLAB_API_VERSION "")
    if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.14.0")
        set(CALS_MATLAB_API_VERSION "R2018a")
    endif ()

    matlab_add_mex(NAME cp_cals_driver
            SRC ${CMAKE_SOURCE_DIR}/matlab/matlab_cp_cals.cpp
            ${CALS_MATLAB_API_VERSION})

    target_link_libraries(cp_cals_driver ${CALS_MATLAB_LINK_LIBS})
endif ()

################################
# Experiments and examples
################################
if (NOT HAVE_MATLAB)
    add_subdirectory(src/examples)
    add_subdirectory(src/experiments)
endif ()
