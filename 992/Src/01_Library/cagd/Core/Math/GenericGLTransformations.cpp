//----------------------------------------------------------------------------------
// File:        Core/Math/GenericGLTransformations.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "GenericGLTransformations.h"

#include <new>

using namespace std;

namespace cagd
{
    //(*@\Green{// default constructor, loads the identity matrix}@*)
    GLTransformation::GLTransformation()
    {
        #pragma omp parallel for
        for (GLint i = 0; i < 16; i++)
        {
            _matrix[i] = ((i % 5 == 0) ? 1.0f : 0.0f);
        }
    }

    //(*@\Green{// matrix addition}@*)
    const GLTransformation GLTransformation::operator +(const GLTransformation &rhs) const
    {
        GLTransformation result(*this);

        #pragma omp parallel for
        for (GLint i = 0; i < 16; i++)
        {
            result._matrix[i] += rhs._matrix[i];
        }

        return result;
    }

    //(*@\Green{// matrix subtraction}@*)
    const GLTransformation GLTransformation::operator -(const GLTransformation &rhs) const
    {
        GLTransformation result(*this);

        #pragma omp parallel for
        for (GLint i = 0; i < 16; i++)
        {
            result._matrix[i] -= rhs._matrix[i];
        }

        return result;
    }

    //(*@\Green{// matrix multiplication}@*)
    const GLTransformation GLTransformation::operator *(const GLTransformation &rhs) const
    {
        GLTransformation result;
        result.loadNullMatrix();

        #pragma omp parallel for
        for (GLint i_j = 0; i_j < 16; i_j++)
        {
            GLint i     = i_j % 4;
            GLint j     = i_j / 4;
            GLint index = 4 * j + i;

            for (GLint k = 0; k < 4; k++)
            {
                result._matrix[index] += _matrix[4 * k + i] * rhs._matrix[4 * j + k];
            }
        }

        return result;
    }

    //(*@\Green{// multiplication by scalar from right}@*)
    const GLTransformation GLTransformation::operator *(const GLfloat &rhs) const
    {
        GLTransformation result(*this);

        #pragma omp parallel for
        for (GLint i = 0; i < 16; i++)
        {
            result._matrix[i] *= rhs;
        }

        return result;
    }

    //(*@\Green{// using homogeneous coordinates, transforms/maps a Cartesian coordinate into another one}@*)
    const Cartesian3 GLTransformation::operator *(const Cartesian3 &rhs) const
    {
        GLfloat hIn[] = {(GLfloat)rhs[0], (GLfloat)rhs[1], (GLfloat)rhs[2], 1.0f};
        GLfloat hOut[] = {0.0f, 0.0f, 0.0f, 0.0f};

        #pragma omp parallel for
        for (GLint j_i = 0; j_i < 16; j_i++)
        {
            GLint j = j_i / 4;
            GLint i = j_i % 4;
            GLfloat member = _matrix[4 * i + j] * hIn[i];
            #pragma omp atomic
            hOut[j] += member;
        }

        return Cartesian3(hOut[0] / hOut[3], hOut[1] / hOut[3], hOut[2] / hOut[3]);
    }

    //(*@\Green{// transforms the given homogeneous coordinate to another one}@*)
    const Homogeneous3 GLTransformation::operator *(const Homogeneous3 &rhs) const
    {
        Homogeneous3 result(0.0f, 0.0f, 0.0f, 0.0f);

        #pragma omp parallel for
        for (GLint j_i = 0; j_i < 16; j_i++)
        {
            GLint j = j_i / 4;
            GLint i = j_i % 4;
            GLfloat member = _matrix[4 * i + j] * rhs[i];
            #pragma omp atomic
            result[j] += member;
        }

        return result;
    }

    //(*@\Green{// division by scalar from right}@*)
    const GLTransformation GLTransformation::operator /(const GLfloat &rhs) const
    {
        assert("GLTransformation::operator /(const GLfloat &rhs) : division by zero!" &&
                rhs != 0.0f);

        GLTransformation result(*this);

        #pragma omp parallel for
        for (GLint i = 0; i < 16; i++)
        {
            result._matrix[i] /= rhs;
        }

        return result;
    }

    //(*@\Green{// add to *this}@*)
    GLTransformation& GLTransformation::operator +=(const GLTransformation &rhs)
    {
        #pragma omp parallel for
        for (GLint i = 0; i < 16; i++)
        {
            _matrix[i] += rhs._matrix[i];
        }

        return *this;
    }

    //(*@\Green{// subtract from *this}@*)
    GLTransformation& GLTransformation::operator -= (const GLTransformation &rhs)
    {
        #pragma omp parallel for
        for (GLint i = 0; i < 16; i++)
        {
            _matrix[i] -= rhs._matrix[i];
        }

        return *this;
    }

    //(*@\Green{// multiplicate *this by a constant}@*)
    GLTransformation& GLTransformation::operator *= (const GLfloat &rhs)
    {
        #pragma omp parallel for
        for (GLint i = 0; i < 16; i++)
        {
            _matrix[i] *= rhs;
        }

        return *this;
    }

    //(*@\Green{// divide *this by a constant}@*)
    GLTransformation& GLTransformation::operator /= (const GLfloat &rhs)
    {
        assert("GLTransformation::operator /=(const GLfloat &constant) : division "
               "by zero!" && rhs != 0.0f);

        #pragma omp parallel for
        for (GLint i = 0; i < 16; i++)
        {
            _matrix[i] /= rhs;
        }

        return *this;
    }

    //(*@\Green{// loads the $4 \times 4$ identity matrix}@*)
    GLvoid GLTransformation::loadIdentity()
    {
        #pragma omp parallel for
        for (GLint i = 0; i < 16; i++)
        {
            _matrix[i] = (i % 5) ? 0.0f : 1.0f;
        }
    }

    //(*@\Green{// loads the $4 \times 4$ null matrix}@*)
    GLvoid GLTransformation::loadNullMatrix()
    {
        #pragma omp parallel for
        for (GLint i = 0; i < 16; i++)
        {
            _matrix[i] = 0.0f;
        }
    }

    //(*@\Green{// returns the transpose of the stored matrix}@*)
    GLTransformation GLTransformation::transpose() const
    {
        GLTransformation result;

        #pragma omp parallel for
        for (GLint r_c = 0; r_c < 16; r_c++)
        {
            GLint r = r_c / 4;
            GLint c = r_c % 4;

            result._matrix[c * 4 + r] = _matrix[r * 4 + c];
        }

        return result;
    }

    //(*@\Green{// auxiliar function that evaluates the determinant of a $ 3\times 3$ matrix}@*)
    GLfloat determinat3x3(
            GLfloat a00, GLfloat a01, GLfloat a02,
            GLfloat a10, GLfloat a11, GLfloat a12,
            GLfloat a20, GLfloat a21, GLfloat a22)
    {
        return (a00 * a11 * a22 + a10 * a21 * a02 + a20 * a01 * a12) -
               (a02 * a11 * a20 + a12 * a21 * a00 + a22 * a01 * a10);
    }

    //(*@\Green{// calculates the determinant of the stored matrix}@*)
    GLfloat GLTransformation::determinant() const
    {
        GLfloat
        result = + _matrix[0] * determinat3x3(
                                _matrix[5], _matrix[ 9], _matrix[13],
                                _matrix[6], _matrix[10], _matrix[14],
                                _matrix[7], _matrix[11], _matrix[15])

                 - _matrix[4] * determinat3x3(
                                _matrix[1], _matrix[ 9], _matrix[13],
                                _matrix[2], _matrix[10], _matrix[14],
                                _matrix[3], _matrix[11], _matrix[15])

                 + _matrix[8] * determinat3x3(
                                _matrix[1], _matrix[5], _matrix[13],
                                _matrix[2], _matrix[6], _matrix[14],
                                _matrix[3], _matrix[7], _matrix[15])

                 - _matrix[12] * determinat3x3(
                                _matrix[1], _matrix[5], _matrix[ 9],
                                _matrix[2], _matrix[6], _matrix[10],
                                _matrix[3], _matrix[7], _matrix[11]);

        return result;
    }

    //(*@\Green{// calculates the inverse of the stored matrix}@*)
    //(*@\Green{// if the stored matrix is singular, a 4 x 4 identity matrix will be returned}@*)
    GLTransformation GLTransformation::inverse(bool *invertible) const
    {
        #define SWAP_ROWS(a, b) { GLfloat *_tmp = a; (a)=(b); (b)=_tmp; }
        #define MATRIX(m,r,c) (m)[(c)*4+(r)]

            GLTransformation result;

            GLfloat wtmp[4][8];
            GLfloat m0, m1, m2, m3, s;
            GLfloat *r0, *r1, *r2, *r3;

            r0 = wtmp[0], r1 = wtmp[1], r2 = wtmp[2], r3 = wtmp[3];

            r0[0] = MATRIX(_matrix, 0, 0);
            r0[1] = MATRIX(_matrix, 0, 1);
            r0[2] = MATRIX(_matrix, 0, 2);
            r0[3] = MATRIX(_matrix, 0, 3);
            r0[4] = 1.0, r0[5] = r0[6] = r0[7] = 0.0f;

            r1[0] = MATRIX(_matrix, 1, 0);
            r1[1] = MATRIX(_matrix, 1, 1);
            r1[2] = MATRIX(_matrix, 1, 2);
            r1[3] = MATRIX(_matrix, 1, 3);
            r1[5] = 1.0, r1[4] = r1[6] = r1[7] = 0.0f;

            r2[0] = MATRIX(_matrix, 2, 0);
            r2[1] = MATRIX(_matrix, 2, 1);
            r2[2] = MATRIX(_matrix, 2, 2);
            r2[3] = MATRIX(_matrix, 2, 3);
            r2[6] = 1.0, r2[4] = r2[5] = r2[7] = 0.0f;

            r3[0] = MATRIX(_matrix, 3, 0);
            r3[1] = MATRIX(_matrix, 3, 1);
            r3[2] = MATRIX(_matrix, 3, 2);
            r3[3] = MATRIX(_matrix, 3, 3);
            r3[7] = 1.0, r3[4] = r3[5] = r3[6] = 0.0f;

            //(*@\Green{// choose pivot -- or die}@*)
            if (abs(r3[0]) > abs(r2[0]))
            {
                SWAP_ROWS(r3, r2);
            }

            if (abs(r2[0]) > abs(r1[0]))
            {
                SWAP_ROWS(r2, r1);
            }

            if (abs(r1[0]) > abs(r0[0]))
            {
                SWAP_ROWS(r1, r0);
            }

            if (0.0f == r0[0])
            {
                if (invertible)
                {
                    *invertible = false;
                }
                return result;
            }

            //(*@\Green{// eliminate first variable}@*)
            m1 = r1[0] / r0[0];
            m2 = r2[0] / r0[0];
            m3 = r3[0] / r0[0];

            s      = r0[1];
            r1[1] -= m1 * s;
            r2[1] -= m2 * s;
            r3[1] -= m3 * s;

            s      = r0[2];
            r1[2] -= m1 * s;
            r2[2] -= m2 * s;
            r3[2] -= m3 * s;

            s      = r0[3];
            r1[3] -= m1 * s;
            r2[3] -= m2 * s;
            r3[3] -= m3 * s;

            s = r0[4];
            if (s != 0.0f)
            {
                r1[4] -= m1 * s;
                r2[4] -= m2 * s;
                r3[4] -= m3 * s;
            }

            s = r0[5];
            if (s != 0.0f)
            {
                r1[5] -= m1 * s;
                r2[5] -= m2 * s;
                r3[5] -= m3 * s;
            }

            s = r0[6];
            if (s != 0.0f)
            {
                r1[6] -= m1 * s;
                r2[6] -= m2 * s;
                r3[6] -= m3 * s;
            }
            s = r0[7];
            if (s != 0.0f)
            {
                r1[7] -= m1 * s;
                r2[7] -= m2 * s;
                r3[7] -= m3 * s;
            }

            //(*@\Green{// choose pivot -- or die}@*)
            if (abs(r3[1]) > abs(r2[1]))
            {
                SWAP_ROWS(r3, r2);
            }

            if (abs(r2[1]) > abs(r1[1]))
            {
                SWAP_ROWS(r2, r1);
            }

            if (0.0f == r1[1])
            {
                if (invertible)
                {
                    *invertible = false;
                }
                return result;
            }

            //(*@\Green{// eliminate second variable}@*)
            m2 = r2[1] / r1[1];
            m3 = r3[1] / r1[1];

            r2[2] -= m2 * r1[2];
            r3[2] -= m3 * r1[2];
            r2[3] -= m2 * r1[3];
            r3[3] -= m3 * r1[3];

            s = r1[4];
            if (0.0f != s)
            {
                r2[4] -= m2 * s;
                r3[4] -= m3 * s;
            }

            s = r1[5];
            if (0.0f != s)
            {
                r2[5] -= m2 * s;
                r3[5] -= m3 * s;
            }

            s = r1[6];
            if (0.0f != s)
            {
                r2[6] -= m2 * s;
                r3[6] -= m3 * s;
            }

            s = r1[7];
            if (0.0f != s)
            {
                r2[7] -= m2 * s;
                r3[7] -= m3 * s;
            }

            //(*@\Green{// choose pivot -- or die}@*)
            if (abs(r3[2]) > abs(r2[2]))
            {
                SWAP_ROWS(r3, r2);
            }

            if (0.0f == r2[2])
            {
                if (invertible)
                {
                    *invertible = false;
                }
                return result;
            }

            //(*@\Green{// eliminate third variable}@*)
            m3 = r3[2] / r2[2];
            r3[3] -= m3 * r2[3], r3[4] -= m3 * r2[4],
            r3[5] -= m3 * r2[5], r3[6] -= m3 * r2[6], r3[7] -= m3 * r2[7];

            //(*@\Green{// last check}@*)
            if (0.0f == r3[3])
            {
                if (invertible)
                {
                    *invertible = false;
                }
                return result;
            }

            s      = 1.0f / r3[3];		  //(*@\Green{// now back substitute row 3}@*)
            r3[4] *= s;
            r3[5] *= s;
            r3[6] *= s;
            r3[7] *= s;

            m2    = r2[3];			  //(*@\Green{// now back substitute row 2}@*)
            s     = 1.0f / r2[2];
            r2[4] = s * (r2[4] - r3[4] * m2), r2[5] = s * (r2[5] - r3[5] * m2),
            r2[6] = s * (r2[6] - r3[6] * m2), r2[7] = s * (r2[7] - r3[7] * m2);

            m1     = r1[3];
            r1[4] -= r3[4] * m1, r1[5] -= r3[5] * m1,
            r1[6] -= r3[6] * m1, r1[7] -= r3[7] * m1;

            m0     = r0[3];
            r0[4] -= r3[4] * m0, r0[5] -= r3[5] * m0,
            r0[6] -= r3[6] * m0, r0[7] -= r3[7] * m0;

            m1    = r1[2];			  //(*@\Green{// now back substitute row 1}@*)
            s     = 1.0f / r1[1];
            r1[4] = s * (r1[4] - r2[4] * m1), r1[5] = s * (r1[5] - r2[5] * m1),
            r1[6] = s * (r1[6] - r2[6] * m1), r1[7] = s * (r1[7] - r2[7] * m1);

            m0     = r0[2];
            r0[4] -= r2[4] * m0, r0[5] -= r2[5] * m0,
            r0[6] -= r2[6] * m0, r0[7] -= r2[7] * m0;

            m0    = r0[1];			  //(*@\Green{// now back substitute row 0}@*)
            s     = 1.0f / r0[0];
            r0[4] = s * (r0[4] - r1[4] * m0), r0[5] = s * (r0[5] - r1[5] * m0),
            r0[6] = s * (r0[6] - r1[6] * m0), r0[7] = s * (r0[7] - r1[7] * m0);

            MATRIX(result._matrix, 0, 0) = r0[4];
            MATRIX(result._matrix, 0, 1) = r0[5];
            MATRIX(result._matrix, 0, 2) = r0[6];
            MATRIX(result._matrix, 0, 3) = r0[7];

            MATRIX(result._matrix, 1, 0) = r1[4];
            MATRIX(result._matrix, 1, 1) = r1[5];
            MATRIX(result._matrix, 1, 2) = r1[6];
            MATRIX(result._matrix, 1, 3) = r1[7];

            MATRIX(result._matrix, 2, 0) = r2[4];
            MATRIX(result._matrix, 2, 1) = r2[5];
            MATRIX(result._matrix, 2, 2) = r2[6];
            MATRIX(result._matrix, 2, 3) = r2[7];

            MATRIX(result._matrix, 3, 0) = r3[4];
            MATRIX(result._matrix, 3, 1) = r3[5];
            MATRIX(result._matrix, 3, 2) = r3[6];
            MATRIX(result._matrix, 3, 3) = r3[7];

            if (invertible)
            {
                *invertible = true;
            }
            return result;

        #undef MATRIX
        #undef SWAP_ROWS
    }

    //(*@\Green{// returns the constant memory address of the stored matrix}@*)
    const GLfloat* GLTransformation::address() const
    {
        return _matrix;
    }

    //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    GLTransformation* GLTransformation::clone() const
    {
        return new (nothrow) GLTransformation(*this);
    }

    //(*@\Green{// virtual default destructor}@*)
    GLTransformation::~GLTransformation() = default;

    //(*@\Green{// multiplicate the transformation matrix rhs by a constant from left}@*)
    const GLTransformation operator *(GLfloat lhs, const GLTransformation &rhs)
    {
        GLTransformation result(rhs);

        #pragma omp parallel for
        for (GLint i = 0; i < 16; i++)
        {
            result._matrix[i] *= lhs;
        }

        return result;
    }

    //(*@\Green{// output to stream}@*)
    std::ostream& operator <<(std::ostream &lhs, const GLTransformation &rhs)
    {
        for (GLint i = 0; i < 16; i++)
        {
            lhs << rhs[i] << " ";
        }

        return lhs;
    }

    //(*@\Green{// input from stream}@*)
    std::istream& operator >>(std::istream &lhs, GLTransformation &rhs)
    {
        for (GLint i = 0; i < 16; i++)
        {
            lhs >> rhs[i];
        }

        return lhs;
    }
}
