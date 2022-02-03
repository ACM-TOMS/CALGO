//----------------------------------------------------------------------------------
// File:        Core/Geometry/Coordinates/Homogeneous3.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef HOMOGENEOUS3_H
#define HOMOGENEOUS3_H

#include <GL/glew.h>

#include "Cartesians3.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <new>

namespace cagd
{
    class Homogeneous3
    {
    private:
        //(*@\Green{// coordinates}@*)
        GLfloat _coord[4];

    public:
        //(*@\Green{// default constructor}@*)
        Homogeneous3();

        //(*@\Green{// special constructor}@*)
        Homogeneous3(GLfloat x, GLfloat y, GLfloat z = 0.0, GLfloat w = 1.0);

        //(*@\Green{// special constructor}@*)
        explicit Homogeneous3(const Cartesian3 &c);

        //(*@\Green{// get components by constant references}@*)
        const GLfloat& operator [](const GLint &rhs) const;
        const GLfloat& x() const;
        const GLfloat& y() const;
        const GLfloat& z() const;
        const GLfloat& w() const;

        //(*@\Green{// get components by non-constant references}@*)
        GLfloat& operator [](const GLint &rhs);
        GLfloat& x();
        GLfloat& y();
        GLfloat& z();
        GLfloat& w();

        //(*@\Green{// sign changing unary operators}@*)
        const Homogeneous3 operator +() const;
        const Homogeneous3 operator -() const;

        //(*@\Green{// addition}@*)
        const Homogeneous3 operator +(const Homogeneous3& rhs) const;

        //(*@\Green{// adds to *this}@*)
        Homogeneous3& operator +=(const Homogeneous3& rhs);

        //(*@\Green{// subtraction}@*)
        const Homogeneous3 operator -(const Homogeneous3& rhs) const;

        //(*@\Green{// subtracts from *this}@*)
        Homogeneous3& operator -=(const Homogeneous3& rhs);

        //(*@\Green{// cross product}@*)
        const Homogeneous3 operator ^(const Homogeneous3& rhs) const;

        //(*@\Green{// cross product, result is stored by *this}@*)
        Homogeneous3& operator ^=(const Homogeneous3& rhs);

        //(*@\Green{// dot product}@*)
        GLfloat operator *(const Homogeneous3& rhs) const;

        //(*@\Green{// scale from right}@*)
        const Homogeneous3 operator *(GLfloat rhs) const;
        const Homogeneous3 operator /(GLfloat rhs) const;

        //(*@\Green{// scale *this}@*)
        Homogeneous3& operator *=(GLfloat rhs);
        Homogeneous3& operator /=(GLfloat rhs);

        //(*@\Green{// returns the Euclidean norm of the represented Cartesian coordinate}@*)
        GLfloat length() const;

        //(*@\Green{// normalizes the represented Cartesian coordinate}@*)
        Homogeneous3& normalize();

        //(*@\Green{// get constant pointer to constant data}@*)
        const GLfloat * address() const;

        //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        Homogeneous3* clone() const;
    };

    //(*@\Green{// default constructor}@*)
    inline Homogeneous3::Homogeneous3()
    {
        _coord[0] = _coord[1] = _coord[2] = 0.0;
        _coord[3] = 1.0;
    }

    //(*@\Green{// special constructor}@*)
    inline Homogeneous3::Homogeneous3(GLfloat x, GLfloat y, GLfloat z, GLfloat w)
    {
        _coord[0] = x;
        _coord[1] = y;
        _coord[2] = z;
        _coord[3] = w;
    }

    //(*@\Green{// special constructor}@*)
    inline Homogeneous3::Homogeneous3(const Cartesian3 &c)
    {
        _coord[0] = (GLfloat)c[0];
        _coord[1] = (GLfloat)c[1];
        _coord[2] = (GLfloat)c[2];
        _coord[3] = 1.0f;
    }

    //(*@\Green{// get components by constant references}@*)
    inline const GLfloat& Homogeneous3::operator [](const GLint &rhs) const
    {
        assert("Homogeneous coordinate index is out of bounds!" && (rhs >= 0 && rhs < 4));
        return _coord[rhs];
    }

    inline const GLfloat& Homogeneous3::x() const
    {
        return _coord[0];
    }

    inline const GLfloat& Homogeneous3::y() const
    {
        return _coord[1];
    }

    inline const GLfloat& Homogeneous3::z() const
    {
        return _coord[2];
    }

    inline const GLfloat& Homogeneous3::w() const
    {
        return _coord[3];
    }

    //(*@\Green{// get components by non-constant references}@*)
    inline GLfloat& Homogeneous3::operator [](const GLint &rhs)
    {
        assert("Homogeneous coordinate index is out of bounds!" && (rhs >= 0 && rhs < 4));
        return _coord[rhs];
    }

    inline GLfloat& Homogeneous3::x()
    {
        return _coord[0];
    }

    inline GLfloat& Homogeneous3::y()
    {
        return _coord[1];
    }

    inline GLfloat& Homogeneous3::z()
    {
        return _coord[2];
    }

    inline GLfloat& Homogeneous3::w()
    {
        return _coord[3];
    }

    //(*@\Green{// sign changing unary operators}@*)
    inline const Homogeneous3 Homogeneous3::operator +() const
    {
        return Homogeneous3(_coord[0], _coord[1], _coord[2], _coord[3]);
    }

    inline const Homogeneous3 Homogeneous3::operator -() const
    {
        return Homogeneous3(-_coord[0], -_coord[1], -_coord[2], _coord[3]);
    }

    //(*@\Green{// addition}@*)
    inline const Homogeneous3 Homogeneous3::operator +(const Homogeneous3& rhs) const
    {
        return Homogeneous3(rhs._coord[3] * _coord[0] + _coord[3] * rhs._coord[0],
                            rhs._coord[3] * _coord[1] + _coord[3] * rhs._coord[1],
                            rhs._coord[3] * _coord[2] + _coord[3] * rhs._coord[2],
                            _coord[3] * rhs._coord[3]);
    }

    //(*@\Green{// add to *this}@*)
    inline Homogeneous3& Homogeneous3::operator +=(const Homogeneous3& rhs)
    {
        _coord[0] = rhs._coord[3] * _coord[0] + _coord[3] * rhs._coord[0];
        _coord[1] = rhs._coord[3] * _coord[1] + _coord[3] * rhs._coord[1];
        _coord[2] = rhs._coord[3] * _coord[2] + _coord[3] * rhs._coord[2];
        _coord[3] = _coord[3] * rhs._coord[3];

        return *this;
    }

    //(*@\Green{// subtraction}@*)
    inline const Homogeneous3 Homogeneous3::operator -(const Homogeneous3& rhs) const
    {
        return Homogeneous3(rhs._coord[3] * _coord[0] - _coord[3] * rhs._coord[0],
                            rhs._coord[3] * _coord[1] - _coord[3] * rhs._coord[1],
                            rhs._coord[3] * _coord[2] - _coord[3] * rhs._coord[2],
                            _coord[3] * rhs._coord[3]);
    }

    //(*@\Green{// subtract from *this}@*)
    inline Homogeneous3& Homogeneous3::operator -=(const Homogeneous3& rhs)
    {
        _coord[0] = rhs._coord[3] * _coord[0] - _coord[3] * rhs._coord[0];
        _coord[1] = rhs._coord[3] * _coord[1] - _coord[3] * rhs._coord[1];
        _coord[2] = rhs._coord[3] * _coord[2] - _coord[3] * rhs._coord[2];
        _coord[3] = _coord[3] * rhs._coord[3];

        return *this;
    }

    //(*@\Green{// cross product}@*)
    inline const Homogeneous3 Homogeneous3::operator ^(const Homogeneous3& rhs) const
    {
        return Homogeneous3(_coord[1] * rhs._coord[2] - _coord[2] * rhs._coord[1],
                            _coord[2] * rhs._coord[0] - _coord[0] * rhs._coord[2],
                            _coord[0] * rhs._coord[1] - _coord[1] * rhs._coord[0],
                            _coord[3] * rhs._coord[3]);
    }

    //(*@\Green{// cross product, result is stored by *this}@*)
    inline Homogeneous3& Homogeneous3::operator ^=(const Homogeneous3& rhs)
    {
        GLfloat x = _coord[1] * rhs._coord[2] - _coord[2] * rhs._coord[1],
                y = _coord[2] * rhs._coord[0] - _coord[0] * rhs._coord[2],
                z = _coord[0] * rhs._coord[1] - _coord[1] * rhs._coord[0],
                w = _coord[3] * rhs._coord[3];

        _coord[0] = x;
        _coord[1] = y;
        _coord[2] = z;
        _coord[3] = w;

        return *this;
    }

    //(*@\Green{// dot product}@*)
    inline GLfloat Homogeneous3::operator *(const Homogeneous3& rhs) const
    {
        if (_coord[3] == 0.0f || rhs._coord[3] == 0.0f)
        {
            return std::numeric_limits<GLfloat>::max();
        }

        return (_coord[0] * rhs._coord[0] +
                _coord[1] * rhs._coord[1] +
                _coord[2] * rhs._coord[2]) / _coord[3] / rhs._coord[3];
    }

    //(*@\Green{// scale from right}@*)
    inline const Homogeneous3 Homogeneous3::operator *(GLfloat rhs) const
    {
        return Homogeneous3(_coord[0] * rhs, _coord[1] * rhs, _coord[2] * rhs, _coord[3]);
    }

    //(*@\Green{// scale from left}@*)
    inline const Homogeneous3 operator *(GLfloat lhs, const Homogeneous3& rhs)
    {
        return Homogeneous3(lhs * rhs[0], lhs * rhs[1], lhs * rhs[2], rhs[3]);
    }

    //(*@\Green{// scale from right}@*)
    inline const Homogeneous3 Homogeneous3::operator /(GLfloat rhs) const
    {
        return Homogeneous3(_coord[0], _coord[1], _coord[2], _coord[3] * rhs);
    }

    //(*@\Green{// scale *this}@*)
    inline Homogeneous3& Homogeneous3::operator *=(GLfloat rhs)
    {
        _coord[0] *= rhs;
        _coord[1] *= rhs;
        _coord[2] *= rhs;

        return *this;
    }

    inline Homogeneous3& Homogeneous3::operator /=(GLfloat rhs)
    {
        _coord[3] *= rhs;

        return *this;
    }

    //(*@\Green{// returns the Euclidean norm of the represented Cartesian coordinate}@*)
    inline GLfloat Homogeneous3::length() const
    {
        return std::sqrt((*this) * (*this));
    }

    //(*@\Green{// normalize}@*)
    inline Homogeneous3& Homogeneous3::normalize()
    {
        GLfloat l = length();

        if (l && l != 1.0)
        {
            *this /= l;
        }

        return *this;
    }

    //(*@\Green{// get constant pointer to constant data}@*)
    inline const GLfloat* Homogeneous3::address() const
    {
        return _coord;
    }

    //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    inline Homogeneous3* Homogeneous3::clone() const
    {
        return new (std::nothrow) Homogeneous3(_coord[0], _coord[1], _coord[2],
                                               _coord[3]);
    }

    //(*@\Green{// output to stream}@*)
    inline std::ostream& operator <<(std::ostream& lhs, const Homogeneous3& rhs)
    {
        return lhs << rhs[0] << " " << rhs[1] << " " << rhs[2] << " " << rhs[3];
    }

    //(*@\Green{// input from stream}@*)
    inline std::istream& operator >>(std::istream& lhs, Homogeneous3& rhs)
    {
        return lhs >> rhs[0] >> rhs[1] >> rhs[2] >> rhs[3];
    }
}

#endif //(*@\Green{// HOMOGENEOUS3\_H}@*)
