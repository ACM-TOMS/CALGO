
SET(ACML_USE_SHARED_LIBRARIES OFF CACHE BOOL "Use the ACLM shared library." )
SET(ACML_USE_OPENMP_LIBRARIES ON CACHE BOOL "Use the OpenMP implementation.")



MACRO(RESET_VAR_IF_CHANGED toChange toCheck)
	IF( NOT ACML_${toCheck}_PREV_VALUE )
		IF( ${toCheck} )
			SET(NEED_TO_CLEAR TRUE)
		ENDIF( ${toCheck} )
	ENDIF( NOT ACML_${toCheck}_PREV_VALUE )


	IF( NOT ${toCheck} )
		IF( ACML_${toCheck}_PREV_VALUE )
			SET(NEED_TO_CLEAR TRUE)
		ENDIF( ACML_${toCheck}_PREV_VALUE )
	ENDIF( NOT ${toCheck} )

	IF( NEED_TO_CLEAR )
		SET(${toChange} "NOTFOUND" CACHE PATH "" FORCE)
	ENDIF( NEED_TO_CLEAR )

	SET(ACML_${toCheck}_PREV_VALUE ${${toCheck}} CACHE INTERNAL "" FORCE)
ENDMACRO(RESET_VAR_IF_CHANGED toChange)


RESET_VAR_IF_CHANGED( ACML_INCLUDE_PATH ACML_USE_OPENMP_LIBRARIES)
RESET_VAR_IF_CHANGED( ACML_INCLUDE_PATH ACML_USE_SHARED_LIBRARIES)
RESET_VAR_IF_CHANGED( ACML ACML_USE_OPENMP_LIBRARIES)
RESET_VAR_IF_CHANGED( ACML ACML_USE_SHARED_LIBRARIES)
RESET_VAR_IF_CHANGED( ACML_SEARCH_PATHS ACML_USE_OPENMP_LIBRARIES)
RESET_VAR_IF_CHANGED( ACML_SEARCH_PATHS ACML_USE_SHARED_LIBRARIES)

SET(ACML_BASE_SEARCH_PATHS "/opt/acml" "C:\\AMD\\acml" "C:\\Program Files\\AMD\\acml")
SET(ACML_SEARCH_PATHS "")

# The various supported versions.  When new versions come out this will need
# to be updated to allow automatic detection.
SET(ACML_VERSIONS 4.2.0 4.1.0 4.0.0 3.5.0 3.6.0 3.6.1)

SET(MP_COMPILER_VERSIONS gfortran64_mp gfortran64_mp_int64 pgi32_mp ifort32_mp)
SET(COMPILER_VERSIONS gfortran64 gfortran64_int64 pgi32 ifort32)

SET(COMPILER_VERSION_TO_USE "a")

IF( ACML_USE_OPENMP_LIBRARIES )
    SET(COMPILER_VERSION_TO_USE ${MP_COMPILER_VERSIONS})
ELSE( ACML_USE_OPENMP_LIBRARIES )
    SET(COMPILER_VERSION_TO_USE ${COMPILER_VERSIONS})
ENDIF( ACML_USE_OPENMP_LIBRARIES )

FOREACH(path_iter ${ACML_BASE_SEARCH_PATHS})
    FOREACH(version_iter ${ACML_VERSIONS})
        FOREACH(compiler_iter ${COMPILER_VERSION_TO_USE})
            LIST(APPEND ACML_SEARCH_PATHS ${path_iter}${version_iter}/${compiler_iter}/include)
        ENDFOREACH(compiler_iter ${COMPILER_VERSION_TO_USE})
    ENDFOREACH(version_iter ${ACML_VERSIONS})	     
ENDFOREACH(path_iter ${ACML_BASE_SEARCH_PATHS})

#MESSAGE(${ACML_SEARCH_PATHS})
FIND_PATH(ACML_INCLUDE_PATH acml.h ${ACML_SEARCH_PATHS} )
SET(ACML_LIB_PATH ${ACML_INCLUDE_PATH}/../lib)
#MESSAGE(${ACML_LIB_PATH})

IF( ACML_USE_SHARED_LIBRARIES )
	IF( ACML_USE_OPENMP_LIBRARIES )
		SET(ACML_LIB_NAMES acml_mp acml_mp_dll libacml_mp_dll)
	ELSE( ACML_USE_OPENMP_LIBRARIES )
		SET(ACML_LIB_NAMES acml acml_dll libacml_dll)
	ENDIF( ACML_USE_OPENMP_LIBRARIES )
ELSE( ACML_USE_SHARED_LIBRARIES )
	IF( ACML_USE_OPENMP_LIBRARIES )
		SET(ACML_LIB_NAMES acml_mp libacml_mp)
	ELSE( ACML_USE_OPENMP_LIBRARIES )
		SET(ACML_LIB_NAMES acml libacml)
	ENDIF( ACML_USE_OPENMP_LIBRARIES )
ENDIF( ACML_USE_SHARED_LIBRARIES )

FIND_LIBRARY( ACML NAMES ${ACML_LIB_NAMES} PATHS ${ACML_LIB_PATH} )

IF( ${CMAKE_GENERATOR} STREQUAL "Unix Makefiles" )
	SET( ACML_TARGET_LINK_LIBRARIES ${ACML} gfortran )
ELSE( ${CMAKE_GENERATOR} STREQUAL "Unix Makefiles" )
	SET( ACML_TARGET_LINK_LIBRARIES ${ACML} )
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Unix Makefiles" )

IF (ACML_INCLUDE_PATH)
  SET(ACML_FOUND ON)
ENDIF (ACML_INCLUDE_PATH)

IF (ACML_FOUND)
  IF (NOT ACML_FIND_QUIETLY)
     MESSAGE(STATUS "Found ACML: ${ACML_INCLUDE_PATH}")
  ENDIF (NOT ACML_FIND_QUIETLY)
ELSE(ACML_FOUND)
  IF (ACML_FIND_REQUIRED)
     MESSAGE(FATAL_ERROR "Could not find ACML")
  ENDIF (ACML_FIND_REQUIRED)
ENDIF (ACML_FOUND)


