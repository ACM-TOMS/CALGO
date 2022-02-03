//----------------------------------------------------------------------------------
// File:        Core/Math/GenericGLTransformations.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef GENERICGLTRANSFORMATIONS_H
#define GENERICGLTRANSFORMATIONS_H

#include "../Geometry/Coordinates/Cartesians3.h"
#include "../Geometry/Coordinates/Homogeneous3.h"

#include <GL/glew.h>
#include <cmath>
#include <iostream>
#include <new>

namespace cagd
{
    class GLTransformation
    {
        //(*@\Green{// multiplicate the transformation matrix rhs by a constant from left}@*)
        friend const GLTransformation operator *(GLfloat lhs,
                                                 const GLTransformation &rhs);

        //(*@\Green{// output to stream}@*)
        friend std::ostream& operator <<(std::ostream &lhs, const GLTransformation &rhs);

        //(*@\Green{// input from stream}@*)
        friend std::istream& operator >>(std::istream &lhs, GLTransformation &rhs);

    protected:
        //(*@\Green{// stores a $4 \times 4$ column-major ordered transformation matrix, i.e.,}@*)
        //(*@\Green{//}@*)
        //(*@\Green{// T = (\_matrix[0] \_matrix[4] \_matrix[ \,8] \_matrix[12]}@*)
        //(*@\Green{// \ \ \ \ \ \ \ \,\_matrix[1] \_matrix[5] \_matrix[ \,9] \_matrix[13]}@*)
        //(*@\Green{// \ \ \ \ \ \ \ \,\_matrix[2] \_matrix[6] \_matrix[10] \_matrix[14]}@*)
        //(*@\Green{// \ \ \ \ \ \ \ \,\_matrix[3] \_matrix[7] \_matrix[11] \_matrix[15])}@*)
        GLfloat _matrix[16];

    public:
        //(*@\Green{// default constructor, loads the identity matrix}@*)
        GLTransformation();

        //(*@\Green{// matrix addition}@*)
        const GLTransformation operator +(const GLTransformation &rhs) const;

        //(*@\Green{// matrix subtraction}@*)
        const GLTransformation operator -(const GLTransformation &rhs) const;

        //(*@\Green{// matrix multiplication}@*)
        const GLTransformation operator *(const GLTransformation &rhs) const;

        //(*@\Green{// multiplication by scalar from right}@*)
        const GLTransformation operator *(const GLfloat &rhs) const;

        //(*@\Green{// using homogeneous coordinates, transforms/maps a Cartesian coordinate into another one}@*)
        const Cartesian3 operator *(const Cartesian3 &rhs) const;

        //(*@\Green{// transforms the given homogeneous coordinate to another one}@*)
        const Homogeneous3 operator *(const Homogeneous3 &rhs) const;

        //(*@\Green{// division by scalar from right}@*)
        const GLTransformation operator /(const GLfloat &rhs) const;

        //(*@\Green{// add to *this}@*)
        GLTransformation& operator +=(const GLTransformation &rhs);

        //(*@\Green{// subtract from *this}@*)
        GLTransformation& operator -=(const GLTransformation &rhs);

        //(*@\Green{// multiplicate *this by a constant}@*)
        GLTransformation& operator *=(const GLfloat &rhs);

        //(*@\Green{// divide *this by a constant}@*)
        GLTransformation& operator /=(const GLfloat &rhs);

        //(*@\Green{// loads the $4 \times 4$ identity matrix}@*)
        GLvoid loadIdentity();

        //(*@\Green{// loads the $4 \times 4$ null matrix}@*)
        GLvoid loadNullMatrix();

        //(*@\Green{// get element by constant reference}@*)
        const GLfloat& operator [](GLint i) const
        {
            assert("The given index is out of bounds!" && (i >= 0 && i < 16));
            return _matrix[i];
        }

        const GLfloat& operator ()(GLint row, GLint column) const
        {
            assert("The given row index is out of bounds!" && (row >= 0 && row < 4));
            assert("The given column index is out of bounds!" &&
                   (column >= 0 && column < 4));
            return _matrix[column * 4 + row];
        }

        //(*@\Green{// get element by non-constant reference}@*)
        GLfloat& operator [](GLint i)
        {
            assert("The given index is out of bounds!" && (i >= 0 && i < 16));
            return _matrix[i];
        }

        GLfloat& operator ()(GLint row, GLint column)
        {
            assert("The given row index is out of bounds!" && (row >= 0 && row < 4));
            assert("The given column index is out of bounds!" &&
                   (column >= 0 && column < 4));
            return _matrix[column * 4 + row];
        }

        //(*@\Green{// returns the transpose of the stored matrix}@*)
        GLTransformation transpose() const;

        //(*@\Green{// calculates the determinant of the stored matrix}@*)
        GLfloat determinant() const;

        //(*@\Green{// calculates the inverse of the stored matrix}@*)
        //(*@\Green{// if the stored matrix is singular, a 4 x 4 identity matrix will be returned}@*)
        GLTransformation inverse(bool *invertible = nullptr) const;

        //(*@\Green{// returns the constant memory address of the stored matrix}@*)
        const GLfloat* address() const; //(*@\label{src:GLTransformation:address:declaration}@*)

        //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        virtual GLTransformation* clone() const;

        //(*@\Green{// virtual default destructor}@*)
        virtual ~GLTransformation();
    };
}

#endif //(*@\Green{// GENERICGLTRANSFORMATIONS\_H}@*)
