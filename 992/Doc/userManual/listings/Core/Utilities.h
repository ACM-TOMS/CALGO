//----------------------------------------------------------------------------------
// File:        Core/Utilities.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef UTILITIES_H
#define UTILITIES_H

#include "Geometry/Coordinates/Colors4.h"

#include <sstream>

namespace cagd
{
    bool        platformIsSupported();
    Color4      coldToHotColormap(GLfloat value, GLfloat min_value, GLfloat max_value); (*@\label{src:Utilities:coldToHotColormap:declaration}@*)
    const char* openGLTypeToString(GLenum type);

    template <typename T>
    std::string toString(const T &value)
    {
        std::ostringstream stream;
        stream << value;
        return stream.str();
    }
}

#endif (*@\Green{// UTILITIES\_H}@*)
