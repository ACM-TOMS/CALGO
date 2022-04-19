
add_library(ExpTempObjects STATIC ExpressionTemplateObjects.cpp)

INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS}
 )


LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})



subdirs(WithoutExpressionTemplates WithExpressionTemplates)
