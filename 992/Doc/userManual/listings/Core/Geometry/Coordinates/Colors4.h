//----------------------------------------------------------------------------------
// File:        Core/Geometry/Coordinates/Colors4.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef COLORS4_H
#define COLORS4_H

#include <GL/glew.h>

#include <algorithm>
#include <cassert>
#include <new>

namespace cagd
{
    class Color4
    {
    private:
        (*@\Green{// red, green, blue and alpha color components}@*)
        GLfloat _data[4];

    public:
        (*@\Green{// default constructor}@*)
        Color4();

        (*@\Green{// special constructor}@*)
        Color4(GLfloat r, GLfloat g, GLfloat b, GLfloat a = 1.0f);

        (*@\Green{// get components by constant references}@*)
        const GLfloat& operator [](const GLint &rhs) const;
        const GLfloat& r() const;
        const GLfloat& g() const;
        const GLfloat& b() const;
        const GLfloat& a() const;

        (*@\Green{// get components by non-constant references}@*)
        GLfloat& operator [](const GLint &rhs);
        GLfloat& r();
        GLfloat& g();
        GLfloat& b();
        GLfloat& a();

        (*@\Green{// overloaded binary arithmetical operators}@*)
        const Color4 operator +(const Color4 &rhs) const;
        const Color4 operator -(const Color4 &rhs) const;
        const Color4 operator *(const Color4 &rhs) const;
        const Color4 operator *(const GLfloat &rhs) const;
        const Color4 operator /(const Color4 &rhs) const;
        const Color4 operator /(const GLfloat &rhs) const;

        Color4& operator +=(const Color4 &rhs);
        Color4& operator -=(const Color4 &rhs);
        Color4& operator *=(const Color4 &rhs);
        Color4& operator *=(const GLfloat &rhs);
        Color4& operator /=(const Color4 &rhs);
        Color4& operator /=(const GLfloat &rhs);

        (*@\Green{// get constant pointer to constant data}@*)
        const GLfloat * address() const;

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        Color4* clone() const;
    };

    (*@\Green{// default constructor}@*)
    inline Color4::Color4()
    {
        _data[0] = _data[1] = _data[2] = 0.0f;
        _data[3] = 1.0f;
    }

    (*@\Green{// special constructor}@*)
    inline Color4::Color4(GLfloat r, GLfloat g, GLfloat b, GLfloat a)
    {
        _data[0] = std::max(0.0f, std::min(r, 1.0f));
        _data[1] = std::max(0.0f, std::min(g, 1.0f));
        _data[2] = std::max(0.0f, std::min(b, 1.0f));
        _data[3] = std::max(0.0f, std::min(a, 1.0f));
    }

    (*@\Green{// get components by non-constant references}@*)
    inline GLfloat& Color4::r()
    {
        return _data[0];
    }

    inline GLfloat& Color4::g()
    {
        return _data[1];
    }

    inline GLfloat& Color4::b()
    {
        return _data[2];
    }

    inline GLfloat& Color4::a()
    {
        return _data[3];
    }

    inline GLfloat& Color4::operator [](const GLint &rhs)
    {
        assert("Color component index is out of bounds!" && (rhs >= 0 && rhs < 4));
        return _data[rhs];
    }

    (*@\Green{// get components by constant references}@*)
    inline const GLfloat& Color4::r() const
    {
        return _data[0];
    }

    inline const GLfloat& Color4::g() const
    {
        return _data[1];
    }

    inline const GLfloat& Color4::b() const
    {
        return _data[2];
    }

    inline const GLfloat& Color4::a() const
    {
        return _data[3];
    }

    inline const GLfloat& Color4::operator [](const GLint &rhs) const
    {
        assert("Color component index is out of bounds!" && (rhs >= 0 && rhs < 4));
        return _data[rhs];
    }

    (*@\Green{// overloaded binary arithmetical operators}@*)
    inline const Color4 Color4::operator +(const Color4 &rhs) const
    {
        return Color4(_data[0] + rhs._data[0],
                      _data[1] + rhs._data[1],
                      _data[2] + rhs._data[2],
                      _data[3] + rhs._data[3]);
    }

    inline const Color4 Color4::operator -(const Color4 &rhs) const
    {
        return Color4(_data[0] - rhs._data[0],
                      _data[1] - rhs._data[1],
                      _data[2] - rhs._data[2],
                      _data[3] - rhs._data[3]);
    }

    inline const Color4 Color4::operator *(const Color4 &rhs) const
    {
        return Color4(_data[0] * rhs._data[0],
                      _data[1] * rhs._data[1],
                      _data[2] * rhs._data[2],
                      _data[3] * rhs._data[3]);
    }

    inline const Color4 Color4::operator *(const GLfloat &rhs) const
    {
        return Color4(_data[0] * rhs,
                      _data[1] * rhs,
                      _data[2] * rhs,
                      _data[3] * rhs);
    }

    inline const Color4 operator *(const GLfloat &lhs, const Color4 &rhs)
    {
        return Color4(lhs * rhs[0], lhs * rhs[1], lhs * rhs[2], lhs * rhs[3]);
    }

    inline const Color4 Color4::operator /(const Color4 &rhs) const
    {
        return Color4(_data[0] / rhs._data[0],
                      _data[1] / rhs._data[1],
                      _data[2] / rhs._data[2],
                      _data[3] / rhs._data[3]);
    }

    inline const Color4 Color4::operator /(const GLfloat &rhs) const
    {
        return Color4(_data[0] / rhs,
                      _data[1] / rhs,
                      _data[2] / rhs,
                      _data[3] / rhs);
    }

    inline Color4& Color4::operator +=(const Color4 &rhs)
    {
        for (int i = 0; i < 4; i++)
        {
            _data[i] += rhs._data[i];
            _data[i] =  std::max(0.0f, std::min(_data[i], 1.0f));
        }

        return *this;
    }

    inline Color4& Color4::operator -=(const Color4 &rhs)
    {
        for (int i = 0; i < 4; i++)
        {
            _data[i] -= rhs._data[i];
            _data[i] =  std::max(0.0f, std::min(_data[i], 1.0f));
        }

        return *this;
    }

    inline Color4& Color4::operator *=(const Color4 &rhs)
    {
        for (int i = 0; i < 4; i++)
        {
            _data[i] *= rhs._data[i];
            _data[i] =  std::max(0.0f, std::min(_data[i], 1.0f));
        }

        return *this;
    }

    inline Color4& Color4::operator *=(const GLfloat &rhs)
    {
        for (int i = 0; i < 4; i++)
        {
            _data[i] *= rhs;
            _data[i] =  std::max(0.0f, std::min(_data[i], 1.0f));
        }

        return *this;
    }

    inline Color4& Color4::operator /=(const Color4 &rhs)
    {
        for (int i = 0; i < 4; i++)
        {
            _data[i] /= rhs[i];
            _data[i] =  std::max(0.0f, std::min(_data[i], 1.0f));
        }

        return *this;
    }

    inline Color4& Color4::operator /=(const GLfloat &rhs)
    {
        for (int i = 0; i < 4; i++)
        {
            _data[i] /= rhs;
            _data[i] =  std::max(0.0f, std::min(_data[i], 1.0f));
        }

        return *this;
    }

    (*@\Green{// get constant pointer to constant data}@*)
    inline const GLfloat* Color4::address() const
    {
        return _data;
    }

    (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    inline Color4* Color4::clone() const
    {
        return new (std::nothrow) Color4(*this);
    }

    (*@\Green{// predefined colors\label{src:Color4:colors:start}}@*)
    namespace colors
    {
        const Color4        black(0.000000f, 0.000000f, 0.000000f, 0.500000f);   (*@\Green{// {\color[rgb]{0.000000, 0.000000, 0.000000}\rule{0.2cm}{0.2cm}}}@*)
        const Color4         gray(0.375000f, 0.375000f, 0.375000f, 0.500000f);   (*@\Green{// {\color[rgb]{0.375000, 0.375000, 0.375000}\rule{0.2cm}{0.2cm}}}@*)
        const Color4   light_gray(0.500000f, 0.500000f, 0.500000f, 0.500000f);   (*@\Green{// {\color[rgb]{0.500000, 0.500000, 0.500000}\rule{0.2cm}{0.2cm}}}@*)
        const Color4       silver(0.708232f, 0.708232f, 0.708232f, 0.500000f);   (*@\Green{// {\color[rgb]{0.708232, 0.708232, 0.708232}\rule{0.2cm}{0.2cm}}}@*)
        const Color4        white(1.000000f, 1.000000f, 1.000000f, 0.500000f);   (*@\Green{// {\color[rgb]{1.000000, 1.000000, 1.000000}\rule{0.2cm}{0.2cm}}}@*)
        const Color4     dark_red(0.474937f, 0.080613f, 0.063188f, 0.500000f);   (*@\Green{// {\color[rgb]{0.474937, 0.080613, 0.063188}\rule{0.2cm}{0.2cm}}}@*)
        const Color4          red(0.851563f, 0.144531f, 0.113281f, 0.500000f);   (*@\Green{// {\color[rgb]{0.851563, 0.144531, 0.113281}\rule{0.2cm}{0.2cm}}}@*)
        const Color4    light_red(1.000000f, 0.276371f, 0.244388f, 0.500000f);   (*@\Green{// {\color[rgb]{1.000000, 0.276371, 0.244388}\rule{0.2cm}{0.2cm}}}@*)
        const Color4        green(0.000000f, 0.570313f, 0.246094f, 0.500000f);   (*@\Green{// {\color[rgb]{0.000000, 0.570313, 0.246094}\rule{0.2cm}{0.2cm}}}@*)
        const Color4  light_green(0.659983f, 0.750576f, 0.457252f, 0.500000f);   (*@\Green{// {\color[rgb]{0.659983, 0.750576, 0.457252}\rule{0.2cm}{0.2cm}}}@*)
        const Color4    dark_blue(0.156863f, 0.086275f, 0.435294f, 0.500000f);   (*@\Green{// {\color[rgb]{0.156863, 0.086275, 0.435294}\rule{0.2cm}{0.2cm}}}@*)
        const Color4         blue(0.156863f, 0.086275f, 0.652941f, 0.500000f);   (*@\Green{// {\color[rgb]{0.156863, 0.086275, 0.652941}\rule{0.2cm}{0.2cm}}}@*)
        const Color4    baby_blue(0.400000f, 0.478431f, 0.701961f, 0.500000f);   (*@\Green{// {\color[rgb]{0.400000, 0.478431, 0.701961}\rule{0.2cm}{0.2cm}}}@*)
        const Color4         cyan(0.000000f, 0.576471f, 0.866667f, 0.500000f);   (*@\Green{// {\color[rgb]{0.000000, 0.576471, 0.866667}\rule{0.2cm}{0.2cm}}}@*)
        const Color4   light_blue(0.456214f, 0.728115f, 1.000000f, 0.500000f);   (*@\Green{// {\color[rgb]{0.456214, 0.728115, 1.000000}\rule{0.2cm}{0.2cm}}}@*)
        const Color4     ice_blue(0.458824f, 0.772549f, 0.941176f, 0.500000f);   (*@\Green{// {\color[rgb]{0.458824, 0.772549, 0.941176}\rule{0.2cm}{0.2cm}}}@*)
        const Color4  dark_purple(0.592157f, 0.270588f, 0.470588f, 0.500000f);   (*@\Green{// {\color[rgb]{0.592157, 0.270588, 0.470588}\rule{0.2cm}{0.2cm}}}@*)
        const Color4       purple(0.710582f, 0.324697f, 0.564721f, 0.500000f);   (*@\Green{// {\color[rgb]{0.710582, 0.324697, 0.564721}\rule{0.2cm}{0.2cm}}}@*)
        const Color4 light_purple(0.850000f, 0.700000f, 1.000000f, 0.500000f);   (*@\Green{// {\color[rgb]{0.850000, 0.700000, 1.000000}\rule{0.2cm}{0.2cm}}}@*)
        const Color4  dark_orange(0.902344f, 0.468750f, 0.089844f, 0.500000f);   (*@\Green{// {\color[rgb]{0.902344, 0.468750, 0.089844}\rule{0.2cm}{0.2cm}}}@*)
        const Color4       orange(1.000000f, 0.705730f, 0.448600f, 0.500000f);   (*@\Green{// {\color[rgb]{1.000000, 0.705730, 0.448600}\rule{0.2cm}{0.2cm}}}@*)
    } (*@\label{src:Color4:colors:end}@*)
}

#endif (*@\Green{// COLORS4\_H}@*)
