# Loki
SET(THIRDPARTY_BUILD_LOKI OFF CACHE BOOL
    "Build TinyXML library from ThirdParty.")

IF (THIRDPARTY_BUILD_LOKI)
    IF (NOT EXISTS ${TPSRC}/loki-0.1.3.tar.bz2)
        file(DOWNLOAD ${TPURL}/loki-0.1.3.tar.bz2 ${TPSRC}/loki-0.1.3.tar.bz2)
    ENDIF()
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E tar jxf ${TPSRC}/loki-0.1.3.tar.bz2
        WORKING_DIRECTORY ${TPSRC}
    )   
ENDIF()

INCLUDE (FindLoki)
