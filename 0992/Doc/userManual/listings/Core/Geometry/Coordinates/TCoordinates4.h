//----------------------------------------------------------------------------------
// File:        Core/Geometry/Coordinates/TCoordinates4.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef TCOORDINATES4_H
#define TCOORDINATES4_H

#include <GL/glew.h>

#include <cassert>
#include <iostream>
#include <new>

namespace cagd
{
    class TCoordinate4
    {
    private:
        (*@\Green{// four dimensional (projective) texture coordinates $(s, t, r, q)$}@*)
        GLfloat _data[4];

    public:
        (*@\Green{// default constructor}@*)
        TCoordinate4();

        (*@\Green{// special constructor}@*)
        TCoordinate4(GLfloat s, GLfloat t, GLfloat r = 0.0, GLfloat q = 1.0);

        (*@\Green{// get components by constant references}@*)
        const GLfloat& operator[](const GLint &rhs) const;
        const GLfloat& s() const;
        const GLfloat& t() const;
        const GLfloat& r() const;
        const GLfloat& q() const;

        (*@\Green{// get components by non-constant references}@*)
        GLfloat& operator[](const GLint &rhs);
        GLfloat& s();
        GLfloat& t();
        GLfloat& r();
        GLfloat& q();

        (*@\Green{// get constant pointer to constant data}@*)
        const GLfloat* address() const;

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        TCoordinate4* clone() const;
    };

    (*@\Green{// default constructor}@*)
    inline TCoordinate4::TCoordinate4()
    {
        _data[0] = _data[1] = _data[2] = 0.0;
        _data[3] = 1.0;
    }

    (*@\Green{// special constructor}@*)
    inline TCoordinate4::TCoordinate4(GLfloat s, GLfloat t, GLfloat r, GLfloat q)
    {
        _data[0] = s;
        _data[1] = t;
        _data[2] = r;
        _data[3] = q;
    }

    (*@\Green{// get components by constant references}@*)
    inline const GLfloat& TCoordinate4::operator [](const GLint &rhs) const
    {
        assert("Texture coordinate index is out of bounds!" && (rhs >= 0 && rhs < 4));
        return _data[rhs];
    }

    inline const GLfloat& TCoordinate4::s() const
    {
        return _data[0];
    }

    inline const GLfloat& TCoordinate4::t() const
    {
        return _data[1];
    }

    inline const GLfloat& TCoordinate4::r() const
    {
        return _data[2];
    }

    inline const GLfloat& TCoordinate4::q() const
    {
        return _data[3];
    }

    (*@\Green{// get components by non-constant references}@*)
    inline GLfloat& TCoordinate4::operator [](const GLint &rhs)
    {
        assert("Texture coordinate index is out of bounds!" && (rhs >= 0 && rhs < 4));
        return _data[rhs];
    }

    inline GLfloat& TCoordinate4::s()
    {
        return _data[0];
    }

    inline GLfloat& TCoordinate4::t()
    {
        return _data[1];
    }

    inline GLfloat& TCoordinate4::r()
    {
        return _data[2];
    }

    inline GLfloat& TCoordinate4::q()
    {
        return _data[3];
    }

    (*@\Green{// get constant pointer to constant data}@*)
    inline const GLfloat* TCoordinate4::address() const
    {
        return _data;
    }

    (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    inline TCoordinate4* TCoordinate4::clone() const
    {
        return new (std::nothrow) TCoordinate4(*this);
    }

    (*@\Green{// overloaded output to stream operator}@*)
    inline std::ostream& operator <<(std::ostream& lhs, const TCoordinate4 &rhs)
    {
        return lhs << rhs.s() << " " << rhs.t() << rhs.r() << " " << rhs.q();
    }

    (*@\Green{// overloaded input from stream operator}@*)
    inline std::istream& operator >>(std::istream& lhs, TCoordinate4 &rhs)
    {
        return lhs >> rhs[0] >> rhs[1] >> rhs[2] >> rhs[3];
    }
}

#endif (*@\Green{// TCOORDINATES4\_H}@*)
