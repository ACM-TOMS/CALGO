//----------------------------------------------------------------------------------
// File:        Core/Geometry/Coordinates/Cartesians3.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef CARTESIANS3_H
#define CARTESIANS3_H

#include <GL/glew.h>

#include <cassert>
#include <cmath>
#include <iostream>

namespace cagd
{
    class Cartesian3
    {
    private:
        (*@\Green{// coordinates}@*)
        GLdouble _coord[3];

    public:
        (*@\Green{// default constructor}@*)
        Cartesian3();

        (*@\Green{// special constructor}@*)
        Cartesian3(const GLdouble &x, const GLdouble &y, const GLdouble &z = 0.0);

        (*@\Green{// get components by constant references}@*)
        const GLdouble& operator [](const GLint &rhs) const;
        const GLdouble& x() const;
        const GLdouble& y() const;
        const GLdouble& z() const;

        (*@\Green{// get components by non-constant references}@*)
        GLdouble& operator [](const GLint &rhs);
        GLdouble& x();
        GLdouble& y();
        GLdouble& z();

        (*@\Green{// sign changing unary operators}@*)
        const Cartesian3 operator +() const;
        const Cartesian3 operator -() const;
		
        (*@\Green{// addition}@*)
        const Cartesian3 operator +(const Cartesian3 &rhs) const;

        (*@\Green{// add to *this}@*)
        Cartesian3& operator +=(const Cartesian3 &rhs);

        (*@\Green{// subtraction}@*)
        const Cartesian3 operator -(const Cartesian3 &rhs) const;

        (*@\Green{// subtract from *this}@*)
        Cartesian3& operator -=(const Cartesian3 &rhs);

        (*@\Green{// cross product $\mathbf{v}_1\left(x_1,y_1,z_1\right)\times\mathbf{v}_2\left(x_2,y_2,z_2\right) = \det\left[\begin{array}{ccc}\mathbf{i}&\mathbf{j}&\mathbf{k}\\x_1&y_1&z_1\\x_2&y_2&z_2\end{array}\right]$}@*)
        const Cartesian3 operator ^(const Cartesian3 &rhs) const;

        (*@\Green{// cross product, result is stored by *this}@*)
        Cartesian3& operator ^=(const Cartesian3 &rhs);

        (*@\Green{// dot product}@*)
        GLdouble operator *(const Cartesian3 &rhs) const;

        (*@\Green{// scale from right}@*)
        const Cartesian3 operator *(const GLdouble &rhs) const;
        const Cartesian3 operator /(const GLdouble &rhs) const;

        (*@\Green{// scale *this}@*)
        Cartesian3& operator *=(const GLdouble &rhs);
        Cartesian3& operator /=(const GLdouble &rhs);

        (*@\Green{// returns the Euclidean norm of the stored vector}@*)
        GLdouble length() const;

        (*@\Green{// normalizes the stored vector}@*)
        Cartesian3& normalize();

        (*@\Green{// get constant pointer to constant data}@*)
        const GLdouble* address() const;

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        Cartesian3* clone() const;

        (*@\Green{// logical operator}@*)
        GLboolean operator != (const GLdouble &rhs) const;
    };

    (*@\Green{// default constructor}@*)
    inline Cartesian3::Cartesian3()
    {
        _coord[0] = _coord[1] = _coord[2] = 0.0;
    }

    (*@\Green{// special constructor}@*)
    inline Cartesian3::Cartesian3(const GLdouble &x, const GLdouble &y, const GLdouble &z)
    {
        _coord[0] = x;
        _coord[1] = y;
        _coord[2] = z;
    }

    (*@\Green{// get components by constant references}@*)
    inline const GLdouble& Cartesian3::operator [](const GLint &rhs) const
    {
        assert("Cartesian coordinate index is out of bounds!" && (rhs >= 0 && rhs < 3));
        return _coord[rhs];
    }

    inline const GLdouble& Cartesian3::x() const
    {
        return _coord[0];
    }

    inline const GLdouble& Cartesian3::y() const
    {
        return _coord[1];
    }

    inline const GLdouble& Cartesian3::z() const
    {
        return _coord[2];
    }

    (*@\Green{// get components by non-constant references}@*)
    inline GLdouble& Cartesian3::operator [](const GLint &rhs)
    {
        assert("Cartesian coordinate index is out of bounds!" && (rhs >= 0 && rhs < 3));
        return _coord[rhs];
    }

    inline GLdouble& Cartesian3::x()
    {
        return _coord[0];
    }

    inline GLdouble& Cartesian3::y()
    {
        return _coord[1];
    }

    inline GLdouble& Cartesian3::z()
    {
        return _coord[2];
    }

    (*@\Green{// sign changing unary operators}@*)
    inline const Cartesian3 Cartesian3::operator +() const
    {
        return Cartesian3(_coord[0], _coord[1], _coord[2]);
    }

    inline const Cartesian3 Cartesian3::operator -() const
    {
        return Cartesian3(-_coord[0], -_coord[1], -_coord[2]);
    }

    (*@\Green{// addition}@*)
    inline const Cartesian3 Cartesian3::operator +(const Cartesian3 &rhs) const
    {
        return Cartesian3(_coord[0] + rhs._coord[0],
                          _coord[1] + rhs._coord[1],
                          _coord[2] + rhs._coord[2]);
    }

    (*@\Green{// add to *this}@*)
    inline Cartesian3& Cartesian3::operator +=(const Cartesian3 &rhs)
    {
        _coord[0] += rhs._coord[0];
        _coord[1] += rhs._coord[1];
        _coord[2] += rhs._coord[2];

        return *this;
    }

    (*@\Green{// subtraction}@*)
    inline const Cartesian3 Cartesian3::operator -(const Cartesian3 &rhs) const
    {
        return Cartesian3(_coord[0] - rhs._coord[0],
                          _coord[1] - rhs._coord[1],
                          _coord[2] - rhs._coord[2]);
    }

    (*@\Green{// subtract from *this}@*)
    inline Cartesian3& Cartesian3::operator -=(const Cartesian3 &rhs)
    {
        _coord[0] -= rhs._coord[0];
        _coord[1] -= rhs._coord[1];
        _coord[2] -= rhs._coord[2];

        return *this;
    }

    (*@\Green{// cross product}@*)
    inline const Cartesian3 Cartesian3::operator ^(const Cartesian3 &rhs) const
    {
        return Cartesian3(_coord[1] * rhs._coord[2] - _coord[2] * rhs._coord[1],
                          _coord[2] * rhs._coord[0] - _coord[0] * rhs._coord[2],
                          _coord[0] * rhs._coord[1] - _coord[1] * rhs._coord[0]);
    }

    (*@\Green{// cross product, result is stored by *this}@*)
    inline Cartesian3& Cartesian3::operator ^=(const Cartesian3 &rhs)
    {
        GLdouble x = _coord[1] * rhs._coord[2] - _coord[2] * rhs._coord[1],
                 y = _coord[2] * rhs._coord[0] - _coord[0] * rhs._coord[2],
                 z = _coord[0] * rhs._coord[1] - _coord[1] * rhs._coord[0];

        _coord[0] = x;
        _coord[1] = y;
        _coord[2] = z;

        return *this;
    }

    (*@\Green{// dot product}@*)
    inline GLdouble Cartesian3::operator *(const Cartesian3 &rhs) const
    {
        return _coord[0] * rhs._coord[0] +
               _coord[1] * rhs._coord[1] +
               _coord[2] * rhs._coord[2];
    }

    (*@\Green{// scale from right}@*)
    inline const Cartesian3 Cartesian3::operator *(const GLdouble &rhs) const
    {
        return Cartesian3(_coord[0] * rhs, _coord[1] * rhs, _coord[2] * rhs);
    }

    inline const Cartesian3 Cartesian3::operator /(const GLdouble &rhs) const
    {
        return Cartesian3(_coord[0] / rhs, _coord[1] / rhs, _coord[2] / rhs);
    }

    (*@\Green{// scale from left}@*)
    inline const Cartesian3 operator *(const GLdouble& lhs, const Cartesian3 &rhs)
    {
        return Cartesian3(lhs * rhs[0], lhs * rhs[1], lhs * rhs[2]);
    }

    (*@\Green{// scale *this}@*)
    inline Cartesian3& Cartesian3::operator *=(const GLdouble &rhs)
    {
        _coord[0] *= rhs;
        _coord[1] *= rhs;
        _coord[2] *= rhs;

        return *this;
    }

    inline Cartesian3& Cartesian3::operator /=(const GLdouble &rhs)
    {
        _coord[0] /= rhs;
        _coord[1] /= rhs;
        _coord[2] /= rhs;

        return *this;
    }

    (*@\Green{// returns the Euclidean norm of the stored vector}@*)
    inline GLdouble Cartesian3::length() const
    {
        return std::sqrt((*this) * (*this));
    }

    (*@\Green{// normalizes the stored vector}@*)
    inline Cartesian3& Cartesian3::normalize()
    {
        GLdouble l = length();

        if (l && l != 1.0)
        {
            *this /= l;
        }

        return *this;
    }

    (*@\Green{// get constant pointer to constant data}@*)
    inline const GLdouble* Cartesian3::address() const
    {
        return _coord;
    }

    (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    inline Cartesian3* Cartesian3::clone() const
    {
        return new (std::nothrow) Cartesian3(_coord[0], _coord[1], _coord[2]);
    }

    (*@\Green{// logical operator}@*)
    inline GLboolean Cartesian3::operator != (const GLdouble &rhs) const
    {
        return (_coord[0] != rhs || _coord[1] != rhs || _coord[2] != rhs);
    }

    (*@\Green{// output to stream}@*)
    inline std::ostream& operator <<(std::ostream& lhs, const Cartesian3 &rhs)
    {
        return lhs << rhs[0] << " " << rhs[1] << " " << rhs[2];
    }

    (*@\Green{// input from stream}@*)
    inline std::istream& operator >>(std::istream& lhs, Cartesian3 &rhs)
    {
        return lhs >> rhs[0] >> rhs[1] >> rhs[2];
    }
}

#endif (*@\Green{// CARTESIANS3\_H}@*)
