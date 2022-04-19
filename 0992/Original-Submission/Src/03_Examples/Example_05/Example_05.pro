TEMPLATE += app
QT       += gui core opengl
CONFIG   += console

# We assume that the compiler is compatible with the C++ 11 standard.
greaterThan(QT_MAJOR_VERSION, 4){
    CONFIG         += c++11
    QT             += widgets
} else {
    QMAKE_CXXFLAGS += -std=c++0x
}

#-------------------------------------------
# Unix/Linux (but not Macintosh) platform...
#-------------------------------------------

unix: !mac {
    message("Unix/Linux (but not Macintosh) platform...")

    # We assume that the user has already installed:
    #
    # - a compiler that supports at least OpenMP 2.0 and is compatible at least with the C++ 11 standard;
    #
    # - the OpenGL Extension Wrangler (GLEW) library and its include file into the
    #   folders /usr/lib/ and /usr/include/GL, respectively (for more details consider
    #   subsection 1.a of the file "../../01_Library/readme.txt");
    #
    # - both the debug and release builds and the include files of our function library
    #   into the folders "/usr/local/lib" and "/usr/local/include/CAGD", respectively (for more
    #   information consider subsection 2.1.a of the file "../../01_Library/readme.txt");

    # for OpenGL and GLEW installed into /usr/lib/libGLEW.so or /usr/lib/glew.lib
    LIBS        += -lGLEW -lGL

    # the default include-path of the debug/release 32/64-bit variants of our library
    INCLUDEPATH += /usr/local/include/CAGD

    # the default library-path of the debug/release 32/64-bit variants of our library
    # is /usr/local/lib
    CONFIG(release, debug|release): {
        contains(QT_ARCH, i386) {
            message("x86 (i.e., 32-bit) release build")
            LIBS += -lcagd32
        } else {
            message("x64 (i.e., 64-bit) release build")
            LIBS += -lcagd64
        }
    } else: CONFIG(debug, debug|release): {
        contains(QT_ARCH, i386) {
            message("x86 (i.e., 32-bit) debug build")
            LIBS += -lcagd32d
        } else {
            message("x64 (i.e., 64-bit) debug build")
            LIBS += -lcagd64d
        }
    }

    # for activating OpenMP
    QMAKE_CXXFLAGS         += -fopenmp
    QMAKE_LFLAGS           += -fopenmp
    LIBS                   += -fopenmp

    # optimization flag for release mode
    QMAKE_CXXFLAGS_RELEASE *= -O2
}

#----------------------
# Macintosh platform...
#----------------------

mac {
    message("Macintosh platform...")

    # We assume that the user has already installed:
    #
    # - a compiler that supports at least OpenMP 2.0 and is compatible at least with the C++ 11 standard;
    #
    # - the OpenGL Extension Wrangler (GLEW) library by calling the command 'brew install glew' as
    #   a super user (for more details consider subsection 1.b of the file "../../01_Library/readme.txt");
    #
    # - both the debug and release builds and the include files of our function library
    #   into the folders "/usr/local/lib" and "/usr/local/include/CAGD", respectively (for more
    #   information consider subsection 2.1.b of the file "../../01_Library/readme.txt");

    #----------------------------------------
    # IMPORTANT (requires user modification):
    #----------------------------------------
    #
    # - in order to provide correct include- and library paths to the OpenGL Extension Wrangler (GLEW) library,
    #   change the letters x, y, z in the next two lines to the corresponding version numbers of the GLEW
    #   library which was installed by using the command 'brew install glew' (for more information consider
    #   subsection 1.b of the file "../../01_Library/readme.txt").
    INCLUDEPATH += "/usr/local/Cellar/glew/x.y.z/include/"
    LIBS        += -L"/usr/local/Cellar/glew/x.y.z/lib/" -lGLEW

    # on MacOS10+ the OpenGL library has to be added as a framework:
    LIBS        += -framework OpenGL

    # the default include-path of the debug/release 32/64-bit variants of our library
    INCLUDEPATH += /usr/local/include/CAGD

    # the default library-path of the debug/release 32/64-bit variants of our library
    # is /usr/local/lib
    CONFIG(release, debug|release): {
        contains(QT_ARCH, i386) {
            message("x86 (i.e., 32-bit) release build")
            LIBS += -lcagd32
        } else {
            message("x64 (i.e., 64-bit) release build")
            LIBS += -lcagd64
        }
    } else: CONFIG(debug, debug|release): {
        contains(QT_ARCH, i386) {
            message("x86 (i.e., 32-bit) debug build")
            LIBS += -lcagd32d
        } else {
            message("x64 (i.e., 64-bit) debug build")
            LIBS += -lcagd64d
        }
    }

    # for activating OpenMP
    QMAKE_CXXFLAGS         += -fopenmp
    QMAKE_LFLAGS           += -fopenmp
    LIBS                   += -fopenmp

    # optimization flag for release mode
    QMAKE_CXXFLAGS_RELEASE *= -O2
}

#--------------------
# Windows platform...
#--------------------

win32 {
    message("Windows platform...")

    # We assume that the user has already installed:
    #
    # - a compiler that supports at least OpenMP 2.0 and is compatible at least with the C++ 11 standard;
    #
    # - the 32- and 64-bit variants of the OpenGL Extension Wrangler (GLEW) library
    #   into the folders "../../00_Dependencies/Lib/GL/x86/" and "../../00_Dependencies/Lib/GL/x64/",
    #   respectively (for more details consider subsection 1.c of the file "../../01_Library/readme.txt");
    #
    # - the include file (glew.h) of the GLEW library into the folder ".../../00_Dependencies/Include/GL"
    #   (for more information consider subsection 1.c of the file "../../01_Library/readme.txt");
    #
    # - both the debug and release builds of our function library into the folder "../../00_Dependencies/Lib/CAGD/"
    #   (if our library was compiled by means of the Windows-specific make-files provided in the folder
    #   "../../01_Library", then the debug and release builds of our function library had already been copied
    #   to the folder "../../00_Dependencies/Lib/CAGD/" -- for more details consider subsection 2.1.c of the
    #   file "../../01_Library/readme.txt");
    #
    # - the include files of our function library into the folder "../../00_Dependencies/Include/CAGD".
    #   (if our library was compiled by means of the Windows-specific make-files provided in the folder
    #   "../../01_Library", then the include files of our function library had already been copied
    #   to the folder "../Dependencies/Include/CAGD/" -- for more details consider subsection 2.1.c of the
    #   file "../../01_Library/readme.txt").

    INCLUDEPATH += "../../00_Dependencies/Include" "../../00_Dependencies/Include/CAGD"
    DEPENDPATH  += "../../00_Dependencies/Include"

    LIBS += -lopengl32
    CONFIG(release, debug|release): {
        contains(QT_ARCH, i386) {
            message("x86 (i.e., 32-bit) release build")
            LIBS += -L"../../00_Dependencies/Lib/CAGD/" -lcagd32
            LIBS += -L"../../00_Dependencies/Lib/GL/x86/" -lglew32
        } else {
            message("x64 (i.e., 64-bit) release build")
            LIBS += -L"../../00_Dependencies/Lib/CAGD/" -lcagd64
            LIBS += -L"../../00_Dependencies/Lib/GL/x64/" -lglew32
        }
    } else: CONFIG(debug, debug|release): {
        contains(QT_ARCH, i386) {
            message("x86 (i.e., 32-bit) debug build")
            LIBS += -L"../../00_Dependencies/Lib/CAGD/" -lcagd32d
            LIBS += -L"../../00_Dependencies/Lib/GL/x86/" -lglew32
        } else {
            message("x64 (i.e., 64-bit) debug build")
            LIBS += -L"../../00_Dependencies/Lib/CAGD/" -lcagd64d
            LIBS += -L"../../00_Dependencies/Lib/GL/x64" -lglew32
        }
    }

    msvc {
      QMAKE_CXXFLAGS         += -openmp -arch:AVX -D "_CRT_SECURE_NO_WARNINGS"
      QMAKE_CXXFLAGS_RELEASE *= -O2
    }
}

SOURCES += \
    main.cpp \
    GUI/MainWindow.cpp \
    GUI/SideWidget.cpp \
    Spaces/SpecializedECSpaces.cpp \
    GUI/GLWidget.cpp

FORMS += \
    GUI/MainWindow.ui \
    GUI/SideWidget.ui

HEADERS += \
    GUI/MainWindow.h \
    GUI/SideWidget.h \
    GUI/GLWidget.h \
    Spaces/SpecializedECSpaces.h
