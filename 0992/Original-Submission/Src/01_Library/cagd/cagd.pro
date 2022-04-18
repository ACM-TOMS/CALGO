QT -= core gui

TEMPLATE  = lib
CONFIG   += staticlib

# We assume that the compiler is compatible with the C++ 11 standard.
greaterThan(QT_MAJOR_VERSION, 4){
    CONFIG         += c++11
} else {
    QMAKE_CXXFLAGS += -std=c++0x
}

# Based on 32|64-bit variant of the underlying operating system
# and on the debug|realease version that has to be compiled, the
# generated static library will be named cagd{32|64}[d].lib,
# where the letter 'd' denotes the debug version.
TARGET = cagd

CONFIG(release, debug|release): {
    contains(QT_ARCH, i386) {
        message("x86 (i.e., 32-bit) release build")
        TARGET = $$join(TARGET,,,32)
    } else {
        message("x64 (i.e., 64-bit) release build")
        TARGET = $$join(TARGET,,,64)
    }
} else: CONFIG(debug, debug|release): {
    contains(QT_ARCH, i386) {
        message("x86 (i.e., 32-bit) debug build")
        TARGET = $$join(TARGET,,,32d)
    } else {
        message("x64 (i.e., 64-bit) debug build")
        TARGET = $$join(TARGET,,,64d)
    }
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
    #   subsection 1.a of the file "../../01_Library/readme.txt").

    # for OpenGL and GLEW installed into /usr/lib/libGLEW.so or /usr/lib/glew.lib
    LIBS        += -lGLEW -lGL

    # for activating OpenMP
    QMAKE_CXXFLAGS         += -fopenmp
    QMAKE_LFLAGS           += -fopenmp
    LIBS                   += -fopenmp

    # optimization flag for release mode
    QMAKE_CXXFLAGS_RELEASE *= -O2

    # If you are running the Qt Creator as a super user, then you can also copy the generated
    # static library file to the folder "/usr/local/lib", from where later it can be added to
    # other projects, e.g., as:
    #
    #   LIBS += -lcagd64d.lib (in 64bit debug mode);
    #   LIBS += -lcagd64.lib  (in 64bit release mode).
    #
    # Note that, in order to automatically copy the target file to the folder "/usr/local/lib",
    # you also have to add a custom build step in Projects mode on the Build Settings tab:
    #
    # - under Build Steps click on the button Add Build Step;
    # - from the pop-up menu select the Custom Process Step option;
    # - write make and install into the text input boxes Command and Arguments, respectively.
    #
    # After these settings you can uncomment the next two lines:
    # target.path            = /usr/local/lib
    # INSTALLS               += target
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
    #   a super user (for more details consider subsection 1.b of the file "../../01_Library/readme.txt").

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
    # for activating OpenMP
    QMAKE_CXXFLAGS         += -fopenmp
    QMAKE_LFLAGS           += -fopenmp
    LIBS                   += -fopenmp

    # optimization flag for release mode
    QMAKE_CXXFLAGS_RELEASE *= -O2

    # If you are running the Qt Creator as a super user, then you can also copy the generated
    # static library file to the folder "/usr/local/lib", from where later it can be added to
    # other projects, e.g., as:
    #
    #   LIBS += -lcagd64d.lib (in 64bit debug mode);
    #   LIBS += -lcagd64.lib  (in 64bit release mode).
    #
    # Note that, in order to automatically copy the target file to the folder "/usr/local/lib",
    # you also have to add a custom build step in Projects mode on the Build Settings tab:
    #
    # - under Build Steps click on the button Add Build Step;
    # - from the pop-up menu select the Custom Process Step option;
    # - write make and install into the text input boxes Command and Arguments, respectively.
    #
    # After these settings you can uncomment the next two lines:
    # target.path            = /usr/local/lib
    # INSTALLS               += target
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
    #   respectively;
    #
    # - the include file (glew.h) of the GLEW library into the folder "Dependencies/Include/GL".
    #
    # Considering the last two list items:
    #
    # - the link http://glew.sourceforge.net/ also provides the pre-compiled 32- and 64-bit release
    #   binaries glew32.{lib|dll} together with their common include file (glew.h);
    # - the user has to download, extract and copy the previously mentioned files as follows:
    #
    #   glew-x.y.z/bin/Release/Win32/glew32.dll --> ../../00_Dependencies/Lib/GL/x86/glew32.dll
    #   glew-x.y.z/bin/Release/x64/glew32.dll   --> ../../00_Dependencies/Lib/GL/x64/glew32.dll
    #
    #   glew-x.y.z/lib/Release/Win32/glew32.lib --> ../../00_Dependencies/Lib/GL/x86/glew32.lib
    #   glew-x.y.z/lib/Release/x64/glew32.lib   --> ../../00_Dependencies/Lib/GL/x64/glew32.lib
    #
    #   glew-x.y.z/include/GL/glew.h            --> ../../00_Dependencies/Include/GL/glew.h
    #
    # where x.y.z denotes the version number of the downloaded function library GLEW.

    msvc {
        INCLUDEPATH             += "../../00_Dependencies/Include"
        DEPENDPATH              += "../../00_Dependencies/Include"
        QMAKE_CXXFLAGS          += -openmp -arch:AVX -D "_CRT_SECURE_NO_WARNINGS"
        QMAKE_CXXFLAGS_RELEASE  *= -O2

        # OpenGL
        LIBS                    += -lopengl32

        # GLEW
        CONFIG(release, debug|release): {
            contains(QT_ARCH, i386) {
                LIBS += -L"../../00_Dependencies/Lib/GL/x86/" -lglew32
            } else {
                LIBS += -L"../../00_Dependencies/Lib/GL/x64/" -lglew32
            }
        } else: CONFIG(debug, debug|release): {
            contains(QT_ARCH, i386) {
                LIBS += -L"../../00_Dependencies/Lib/GL/x86/" -lglew32
            } else {
                LIBS += -L"../../00_Dependencies/Lib/GL/x64" -lglew32
            }
        }
    }
}

HEADERS += \
    Core/Geometry/Coordinates/Cartesians3.h \
    Core/Geometry/Coordinates/Colors4.h \
    Core/Geometry/Coordinates/Homogeneous3.h \
    Core/Geometry/Coordinates/TCoordinates4.h \
    Core/Geometry/Curves/GenericCurves3.h \
    Core/Geometry/Curves/LinearCombinations3.h \
    Core/Geometry/Surfaces/Lights.h \
    Core/Geometry/Surfaces/Materials.h \
    Core/Geometry/Surfaces/TensorProductSurfaces3.h \
    Core/Geometry/Surfaces/TriangleMeshes3.h \
    Core/Geometry/Surfaces/TriangularFaces.h \
    Core/Math/Constants.h \
    Core/Math/GenericGLTransformations.h \
    Core/Math/Matrices.h \
    Core/Math/PascalTriangles.h \
    Core/Math/RealMatrices.h \
    Core/Math/RealMatrixDecompositions.h \
    Core/Math/SpecialGLTransformations.h \
    Core/Shaders/ShaderPrograms.h \
    Core/SmartPointers/CheckingPolicies.h \
    Core/SmartPointers/ImplicitConversionPolicies.h \
    Core/SmartPointers/OwnershipPolicies.h \
    Core/SmartPointers/SmartPointers.h \
    Core/SmartPointers/SpecializedSmartPointers.h \
    Core/SmartPointers/StaticChecks.h \
    Core/SmartPointers/StoragePolicies.h \
    Core/SmartPointers/TypeSelectors.h \
    Core/Exceptions.h \
    Core/Utilities.h \
    EC/BCurves3.h \
    EC/BSurfaces3.h \
    EC/CharacteristicPolynomials.h \
    EC/ECSpaces.h

SOURCES += \
    Core/Geometry/Curves/GenericCurves3.cpp \
    Core/Geometry/Curves/LinearCombinations3.cpp \
    Core/Geometry/Surfaces/Lights.cpp \
    Core/Geometry/Surfaces/Materials.cpp \
    Core/Geometry/Surfaces/TensorProductSurfaces3.cpp \
    Core/Geometry/Surfaces/TriangleMeshes3.cpp \
    Core/Math/GenericGLTransformations.cpp \
    Core/Math/PascalTriangles.cpp \
    Core/Math/RealMatrices.cpp \
    Core/Math/RealMatrixDecompositions.cpp \
    Core/Math/SpecialGLTransformations.cpp \
    Core/Shaders/ShaderPrograms.cpp \
    Core/Utilities.cpp \
    EC/BCurves3.cpp \
    EC/BSurfaces3.cpp \
    EC/CharacteristicPolynomials.cpp \
    EC/ECSpaces.cpp


