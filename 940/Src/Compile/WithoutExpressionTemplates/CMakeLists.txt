
set(names
  3param_no_transform
  2param
  3param_associative
  3param_commutative
  large_balanced_associative
  loop_fusion
  multipleVectorCommutative
  multiple_associative
)

foreach(i ${names})
  ADD_LIBRARY(${i} STATIC ../${i}.cpp )
  TARGET_LINK_LIBRARIES( ${i}
    ${NEKTAR++_LIBRARIES}
    ${Boost_TIMER_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_CHRONO_LIBRARY}
    ExpTempObjects
  )

  IF( CMAKE_HOST_UNIX )
    TARGET_LINK_LIBRARIES(${i} rt)
  ENDIF()
endforeach()


