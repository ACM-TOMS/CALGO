

SET(Sources
    Armadillo.cpp
    Blitz.cpp
    Eigen.cpp
    Genial.cpp
    HandCoded.cpp
    main.cpp
    NektarMoveConstructors.cpp
    NektarWithExpTemplates.cpp
    NektarWithoutExpTemplates.cpp
    Sets.cpp
    MET.cpp
    UBlas.cpp
)

SET(Headers
    ConjugateGradient.hpp
    ITest.h
    HandCoded.h
    MatrixMultiplication.h
    NektarMoveConstructors.h
    MET.h
    Sets.hpp
    Stat.h
    TimingTestRunner.h
    UBlas.h
    VectorAdditionTests.h
)

#TODO - add -msse and -msse2 to the command line for gcc builds.

SET(Name "ExpressionTemplatePerformance")

INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/ThirdParty/blitz-0.9
    ${CMAKE_SOURCE_DIR}/ThirdParty/met
    ${CMAKE_SOURCE_DIR}/ThirdParty/met/common
    ${CMAKE_SOURCE_DIR}/ThirdParty/freepooma-2.4.1/src
    ${CMAKE_SOURCE_DIR}/ThirdParty/genial-2.1.0/genial
    ${CMAKE_SOURCE_DIR}/ThirdParty/FLENS-2009-06-12-bugfix
    #${CMAKE_SOURCE_DIR}/ThirdParty/armadillo-1.2.0/include
    ${CMAKE_SOURCE_DIR}/ThirdParty/armadillo-3.6.3/include
 )


LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})


ADD_EXECUTABLE(${Name} ${Sources} ${Headers})

TARGET_LINK_LIBRARIES( ${Name}
    ${NEKTAR++_LIBRARIES}
    ${Boost_TIMER_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_CHRONO_LIBRARY}
)

IF( CMAKE_HOST_UNIX )
    TARGET_LINK_LIBRARIES(${Name} rt)
ENDIF()

INSTALL(TARGETS ${Name} 
		RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}) 
