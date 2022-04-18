
SET(Sources
ExpressionTemplateObjects.cpp
main.cpp
TestAliasing.cpp
TestAssociativeTransform.cpp
TestBackwardsInverseTransform.cpp
TestCommutativeTransform.cpp
TestFindACNodesToSwap.cpp
TestForwardInverseTransform.cpp
TestInvertNode.cpp
TestMatrixBlasOptimizations.cpp
TestMatrixBufferSize.cpp
TestMatrixVectorOperations.cpp
TestPushDownUnaryNodes.cpp
TestRemoveUnecessaryTemporaries.cpp
TestTemporaryCount.cpp
TestUnaryOperations.cpp
TestUnrolledVectors.cpp
)

SET(Headers
    ExpressionTemplateObjects.h
)

LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})

SET(Name "ExpressionTemplateTesting")

INCLUDE_DIRECTORIES(../Src ${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS} )

ADD_EXECUTABLE(${Name} ${Sources} ${Headers})

TARGET_LINK_LIBRARIES( ${Name}
	${NEKTAR++_LIBRARIES}
)

INSTALL(TARGETS ${Name} 
		RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}) 