//----------------------------------------------------------------------------------
// File:        Core/Geometry/Surfaces/TriangularFaces.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef TRIANGULARFACES_H
#define TRIANGULARFACES_H

#include <GL/glew.h>

#include <cassert>
#include <iostream>
#include <new>

namespace cagd
{
    class TriangularFace
    {
    private:
        GLint _node[3];

    public:
        //(*@\Green{// default constructor}@*)
        TriangularFace();

        //(*@\Green{// get node identifier by constant reference}@*)
        const GLint& operator [](GLint i) const;

        //(*@\Green{// get node identifier by non-constant reference}@*)
        GLint& operator [](GLint i);

        //(*@\Green{// get constant pointer to constant data}@*)
        const GLint* address() const;

        //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        TriangularFace* clone() const;
    };

    //(*@\Green{// default constructor}@*)
    inline TriangularFace::TriangularFace()
    {
        _node[0] = _node[1] = _node[2] = -1;
    }

    //(*@\Green{// get node identifier by constant reference}@*)
    inline const GLint& TriangularFace::operator [](GLint i) const
    {
        assert("The given node index is out of bounds!" && (i >= 0 && i < 3));
        return _node[i];
    }

    //(*@\Green{// get node identifier by non-constant reference}@*)
    inline GLint& TriangularFace::operator [](GLint i)
    {
        assert("The given node index is out of bounds!" && (i >= 0 && i < 3));
        return _node[i];
    }

    //(*@\Green{// get constant pointer to constant data}@*)
    inline const GLint* TriangularFace::address() const
    {
        return _node;
    }

    //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    inline TriangularFace* TriangularFace::clone() const
    {
        return new (std::nothrow) TriangularFace(*this);
    }

    //(*@\Green{// overloaded output to stream operator}@*)
    inline std::ostream& operator <<(std::ostream& lhs, const TriangularFace& rhs)
    {
        lhs << 3;

        for (GLint i = 0; i < 3; i++)
        {
            lhs  << " " << rhs[i];
        }

        return lhs;
    }

    //(*@\Green{// overloaded input from stream operator}@*)
    inline std::istream& operator >>(std::istream& lhs, TriangularFace& rhs)
    {
        GLint vertex_count;
        return lhs >> vertex_count >> rhs[0] >> rhs[1] >> rhs[2];
    }
}

#endif //(*@\Green{// TRIANGULARFACES\_H}@*)
