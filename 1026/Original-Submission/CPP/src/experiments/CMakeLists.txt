if (WITH_DIAGNOSTICS GREATER_EQUAL 1)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif ()

add_executable(experiments experiments.cpp experiments_utils.cpp)
add_compilation_flags(experiments "${comp_flags}")
target_link_libraries(experiments PUBLIC ${CP_CALS})

add_executable(benchmark_cals_mttkrp benchmark_cals_mttkrp.cpp)
add_compilation_flags(benchmark_cals_mttkrp "${comp_flags}")
target_link_libraries(benchmark_cals_mttkrp PUBLIC ${CP_CALS})

add_executable(evaluate peak_evaluator.cpp)
add_compilation_flags(evaluate "${comp_flags}")
target_link_libraries(evaluate PUBLIC ${CP_CALS})

add_executable(benchmark_other_mttkrp benchmark_other_mttkrp.cpp)
add_compilation_flags(benchmark_other_mttkrp "${comp_flags}")
target_link_libraries(benchmark_other_mttkrp PUBLIC ${CP_CALS})

#############################################################################################
# Letter related experiments
add_executable(experiments_letter experiments_letter.cpp experiments_utils.cpp)
target_link_libraries(experiments_letter PUBLIC ${CP_CALS})

set(WITH_CTF "OFF")
set(WITH_PLANC "OFF")
if (WITH_EXPERIMENTS)
  set(WITH_CTF "ON")
  set(WITH_PLANC "ON")
endif()

if (WITH_CTF)
  set(CTF_PATH "${CMAKE_SOURCE_DIR}/extern/ctf")
    find_package(MPI REQUIRED)

    target_compile_definitions(benchmark_other_mttkrp PUBLIC WITH_CTF=1)
    target_include_directories(benchmark_other_mttkrp PUBLIC "${CTF_PATH}/include")
    target_link_directories(benchmark_other_mttkrp PUBLIC "${CTF_PATH}/lib_shared")
    target_link_libraries(benchmark_other_mttkrp PUBLIC ctf MPI::MPI_CXX)
endif ()

if (WITH_PLANC)
  set(PLANC_PATH "${CMAKE_SOURCE_DIR}/extern/planc")
  set(ARMA_PATH "${CMAKE_SOURCE_DIR}/extern/armadillo")

  include_directories(${PLANC_PATH})
  target_compile_definitions(benchmark_other_mttkrp PUBLIC WITH_PLANC=1)
  target_include_directories(benchmark_other_mttkrp PUBLIC "${PLANC_PATH}/dimtree")
  target_include_directories(benchmark_other_mttkrp PUBLIC "${ARMA_PATH}/include")
endif ()