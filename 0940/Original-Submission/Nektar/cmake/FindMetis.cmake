
SET(METIS_SEARCH_PATHS 
	${CMAKE_SOURCE_DIR}/ThirdParty/modmetis-4.0/
	${CMAKE_SOURCE_DIR}/ThirdParty/modmetis-4.0/build/
	${CMAKE_SOURCE_DIR}/../ThirdParty/modmetis-4.0/
	${CMAKE_SOURCE_DIR}/../ThirdParty/modmetis-4.0/build
    ${CMAKE_SOURCE_DIR}/ThirdParty/dist/lib 
    ${CMAKE_SOURCE_DIR}/../ThirdParty/dist/lib)

FIND_LIBRARY(METIS_LIB NAMES modmetis PATHS ${METIS_SEARCH_PATHS})


SET(METIS_FOUND FALSE)
IF (METIS_LIB)
  SET(METIS_FOUND TRUE)
  MARK_AS_ADVANCED(METIS_LIB)
ENDIF (METIS_LIB)

IF (METIS_FOUND)
  IF (NOT METIS_LIB_FIND_QUIETLY)
     MESSAGE(STATUS "Found Metis")
  ENDIF (NOT METIS_LIB_FIND_QUIETLY)
ELSE(METIS_FOUND)
  IF (METIS_LIB_FIND_REQUIRED)
     MESSAGE(FATAL_ERROR "Could not find Metis")
  ENDIF (METIS_LIB_FIND_REQUIRED)
ENDIF (METIS_FOUND)
