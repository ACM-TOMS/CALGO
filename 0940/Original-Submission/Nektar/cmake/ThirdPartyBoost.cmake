SET(THIRDPARTY_BUILD_BOOST OFF CACHE BOOL
    "Build Boost libraries")

IF (THIRDPARTY_BUILD_BOOST)
    INCLUDE(ExternalProject)
   
    IF (NOT WIN32)
        # Only build the libraries we need
        SET(BOOST_LIB_LIST --with-system --with-iostreams --with-filesystem 
                           --with-program_options --with-date_time --with-thread)
                           
        # We need -fPIC for 64-bit builds
        IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
            SET(BOOST_FLAGS cxxflags=-fPIC cflags=-fPIC linkflags=-fPIC)
        ENDIF ()
        
        # Build Boost
        IF (APPLE)
            SET(TOOLSET darwin)
        ELSE(APPLE)
            SET(TOOLSET gcc)
        ENDIF(APPLE)
        EXTERNALPROJECT_ADD(
            boost
            PREFIX ${TPSRC}
            URL ${TPURL}/boost_1_49_0.tar.bz2
            URL_MD5 "0d202cb811f934282dea64856a175698"
            DOWNLOAD_DIR ${TPSRC}
            CONFIGURE_COMMAND ./bootstrap.sh --prefix=${TPSRC}/dist
            BUILD_COMMAND NO_BZIP2=1 ./b2
                            variant=release
                            link=shared 
                            include=${TPSRC}/dist/include
                            linkflags="-L${TPSRC}/dist/lib"
                            ${BOOST_FLAGS} ${BOOST_LIB_LIST} 
                            --layout=system toolset=${TOOLSET} install
            INSTALL_COMMAND ""
            BUILD_IN_SOURCE 1
        )

        # If building ThirdParty zlib, force zlib build before boost
        IF (THIRDPARTY_BUILD_ZLIB)
            ADD_DEPENDENCIES(boost zlib)
        ENDIF(THIRDPARTY_BUILD_ZLIB)

        # Set up CMake variables
        SET(Boost_DATE_TIME_LIBRARY boost_date_time)
        SET(Boost_DATE_TIME_LIBRARY_DEBUG boost_date_time)
        SET(Boost_DATE_TIME_LIBRARY_RELEASE boost_date_time)
        SET(Boost_FILESYSTEM_LIBRARY boost_filesystem)
        SET(Boost_FILESYSTEM_LIBRARY_DEBUG boost_filesystem)
        SET(Boost_FILESYSTEM_LIBRARY_RELEASE boost_filesystem)
        SET(Boost_IOSTREAMS_LIBRARY boost_iostreams)
        SET(Boost_IOSTREAMS_LIBRARY_DEBUG boost_iostreams)
        SET(Boost_IOSTREAMS_LIBRARY_RELEASE boost_iostreams)
        SET(Boost_PROGRAM_OPTIONS_LIBRARY boost_program_options)
        SET(Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG boost_program_options)
        SET(Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE boost_program_options)
        SET(Boost_SYSTEM_LIBRARY boost_system)
        SET(Boost_SYSTEM_LIBRARY_DEBUG boost_system)
        SET(Boost_SYSTEM_LIBRARY_RELEASE boost_system)
        SET(Boost_THREAD_LIBRARY boost_thread)
        SET(Boost_THREAD_LIBRARY_DEBUG boost_thread)
        SET(Boost_THREAD_LIBRARY_RELEASE boost_thread)
        SET(Boost_INCLUDE_DIRS ${TPSRC}/dist/include 
                               ${TPSRC}/dist/include/boost-1_49)
        SET(Boost_LIBRARY_DIRS ${TPSRC}/dist/lib)
        LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
    ELSE ()
        EXTERNALPROJECT_ADD(
            boost
            PREFIX ${TPSRC}
            URL ${TPURL}/boost_1_49_0.tar.bz2
            URL_MD5 "0d202cb811f934282dea64856a175698"
            DOWNLOAD_DIR ${TPSRC}
            CONFIGURE_COMMAND bootstrap.bat --prefix=${TPSRC}/boost
            BUILD_COMMAND b2 --layout=system install
            INSTALL_COMMAND ""
            BUILD_IN_SOURCE 1
        )
    ENDIF ()
ELSE (THIRDPARTY_BUILD_BOOST)
    SET(Boost_DEBUG 1)
    SET(Boost_USE_MULTITHREAD ON)
    SET(Boost_ADDITIONAL_VERSIONS "1.49" "1.49.0" "1.48" "1.48.0" "1.47.0" "1.47" "1.46" "1.46. 1" "1.40" "1.40.0" "1.35.0" "1.35")
    SET(Boost_NO_BOOST_CMAKE ON)

    IF( NOT BOOST_ROOT )
        #If the user has not set BOOST_ROOT, look in a couple common places first.
        SET(BOOST_ROOT ${CMAKE_SOURCE_DIR}/ThirdParty/boost)
        FIND_PACKAGE( Boost COMPONENTS thread iostreams zlib date_time filesystem system program_options)
        SET(BOOST_ROOT ${CMAKE_SOURCE_DIR}/../ThirdParty/boost)
        FIND_PACKAGE( Boost COMPONENTS thread iostreams zlib date_time filesystem system program_options)
        SET(BOOST_ROOT ${CMAKE_SOURCE_DIR}/ThirdParty/dist)
        FIND_PACKAGE( Boost COMPONENTS thread iostreams zlib date_time filesystem system program_options)
    ELSE()
        FIND_PACKAGE( Boost COMPONENTS thread iostreams zlib date_time filesystem system program_options)
    ENDIF()
ENDIF (THIRDPARTY_BUILD_BOOST)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
