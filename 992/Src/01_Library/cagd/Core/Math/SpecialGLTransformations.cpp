//----------------------------------------------------------------------------------
// File:        Core/Math/SpecialGLTransformations.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "SpecialGLTransformations.h"

namespace cagd
{
    //(*@\Green{// Implementation of translation matrices.}@*)

    //(*@\Green{// default/special constructor}@*)
    Translate::Translate(GLfloat ux, GLfloat uy, GLfloat uz): GLTransformation()
    {
        _matrix[12] = ux;
        _matrix[13] = uy;
        _matrix[14] = uz;
    }

    //(*@\Green{// special constructor, the given direction vector will be normalized}@*)
    Translate::Translate(Cartesian3 &direction, GLfloat distance): GLTransformation()
    {
        direction.normalize();
        _matrix[12] = (GLfloat)direction[0] * distance;
        _matrix[13] = (GLfloat)direction[1] * distance;
        _matrix[14] = (GLfloat)direction[2] * distance;
    }

    //(*@\Green{// setters}@*)
    GLvoid Translate::setXDirectionalUnits(GLfloat ux)
    {
        _matrix[12] = ux;
    }

    GLvoid Translate::setYDirectionalUnits(GLfloat uy)
    {
        _matrix[13] = uy;
    }

    GLvoid Translate::setZDirectionalUnits(GLfloat uz)
    {
        _matrix[14] = uz;
    }

    //(*@\Green{// redefined clone function required by smart pointers based on the deep copy ownership policy}@*)
    Translate* Translate::clone() const
    {
        return new Translate(*this);
    }

    //(*@\Green{// Implementation of scaling transformation.}@*)

    //(*@\Green{// default/special constructor}@*)
    Scale::Scale(GLfloat sx, GLfloat sy, GLfloat sz): GLTransformation()
    {
        _matrix[ 0] = sx;
        _matrix[ 5] = sy;
        _matrix[10] = sz;
    }

    //(*@\Green{// sets the scaling factors along axes $x$, $y$ and $z$}@*)
    GLvoid Scale::setScalingFactors(GLfloat sx, GLfloat sy, GLfloat sz)
    {
        _matrix[ 0] = sx;
        _matrix[ 5] = sy;
        _matrix[10] = sz;
    }

    //(*@\Green{// redefined clone function required by smart pointers based on the deep copy ownership policy}@*)
    Scale* Scale::clone() const
    {
        return new Scale(*this);
    }

    //(*@\Green{// Implementation of rotation matrices.}@*)

    //(*@\Green{// default/special constructor}@*)
    Rotate::Rotate(const Cartesian3 &direction, GLfloat angle_in_radians):
        GLTransformation(),
        _direction(direction), _angle_in_radians(angle_in_radians)
    {
        //(*@\Green{// rotation $R$ will correspond to the given axis of unit direction and the rotation angle given in radians}@*)
        //(*@\Green{//}@*)
        //(*@\Green{// we will use the formula $R = I_3 - \sin(\text{angle\_in\_radians})\cdot S + (1.0 - \cos(\text{angle\_in\_radians})) \cdot S^2$,}@*)
        //(*@\Green{// where $I_3$ is the $3\times 3$ identity matrix, while $S = \left[\begin{array}{rrr}0&\text{direction}[2]&-\text{direction}[1]\\-\text{direction}[2]&0&\text{direction}[0]\\\text{direction}[1]&-\text{direction}[0]&0\end{array}\right]$}@*)
        //(*@\Green{// is a skew-symmetric matrix}@*)
        _direction.normalize();

         GLTransformation S;

         S.loadNullMatrix();
         S[4] =  (GLfloat)_direction[2];
         S[8] = -(GLfloat)_direction[1];
         S[1] = -(GLfloat)_direction[2];
         S[9] =  (GLfloat)_direction[0];
         S[2] =  (GLfloat)_direction[1];
         S[6] = -(GLfloat)_direction[0];

         GLTransformation S2 = S * S;

         S  *= -sin(_angle_in_radians);
         S2 *= 1.0f - cos(_angle_in_radians);

         *this += S;
         *this += S2;
    }

    //(*@\Green{// sets the unit direction vector of the rotation axis}@*)
    GLvoid Rotate::setDirection(const Cartesian3 &direction)
    {
        _direction = direction;
        _direction.normalize();

        GLTransformation S;

        S.loadNullMatrix();
        S[4] =  (GLfloat)_direction[2];
        S[8] = -(GLfloat)_direction[1];
        S[1] = -(GLfloat)_direction[2];
        S[9] =  (GLfloat)_direction[0];
        S[2] =  (GLfloat)_direction[1];
        S[6] = -(GLfloat)_direction[0];

        GLTransformation S2 = S * S;

        S  *= -sin(_angle_in_radians);
        S2 *= 1.0f - cos(_angle_in_radians);

        loadIdentity();

        *this += S;
        *this += S2;
    }

    //(*@\Green{// sets the rotation angle in radians}@*)
    GLvoid Rotate::setAngle(GLfloat angle_in_radians)
    {
        _angle_in_radians = angle_in_radians;

        GLTransformation S;

        S.loadNullMatrix();
        S[4] =  (GLfloat)_direction[2];
        S[8] = -(GLfloat)_direction[1];
        S[1] = -(GLfloat)_direction[2];
        S[9] =  (GLfloat)_direction[0];
        S[2] =  (GLfloat)_direction[1];
        S[6] = -(GLfloat)_direction[0];

        GLTransformation S2 = S * S;

        S  *= -sin(_angle_in_radians);
        S2 *= 1.0f - cos(_angle_in_radians);

        loadIdentity();

        *this += S;
        *this += S2;
    }

    //(*@\Green{// redefined clone function required by smart pointers based on the deep copy ownership policy}@*)
    Rotate* Rotate::clone() const
    {
        return new Rotate(*this);
    }

    //(*@\Green{// Implementation of perspective projection matrices.}@*)

    //(*@\Green{// special constructor}@*)
    PerspectiveProjection::PerspectiveProjection(
            GLfloat aspect, GLfloat fov, GLfloat near, GLfloat far):
        _aspect(aspect), _fov(fov), _f(1.0f / tan(_fov / 2.0f)), _near(near), _far(far)
    {
        _matrix[ 0] =  _f / _aspect;
        _matrix[ 5] =  _f;
        _matrix[10] =  (_near + _far) / (_near - _far);
        _matrix[11] = -1.0f;
        _matrix[14] =  2.0f * _near * _far / (_near - _far);
        _matrix[15] =  0.0f;
    }

    //(*@\Green{// setters}@*)
    GLvoid PerspectiveProjection::setAspectRatio(GLfloat aspect)
    {
        _aspect    = aspect;
        _matrix[0] = _f / _aspect;
    }

    GLvoid PerspectiveProjection::setFieldOfView(GLfloat fov)
    {
        _fov       = fov;
        _f         = 1.0f / tan(_fov / 2.0f);
        _matrix[0] = _f / _aspect;
        _matrix[5] = _f;
    }

    GLvoid PerspectiveProjection::setNearClippingPlaneDistance(GLfloat near)
    {
        _near       = near;
        _matrix[10] = (_near + _far) / (_near - _far);
        _matrix[14] = 2.0f * _near * _far / (_near - _far);
    }

    GLvoid PerspectiveProjection::setFarClippingPlaneDistance(GLfloat far)
    {
        _far        = far;
        _matrix[10] = (_near + _far) / (_near - _far);
        _matrix[14] = 2.0f * _near * _far / (_near - _far);
    }

    //(*@\Green{// redefined clone function required by smart pointers based on the deep copy ownership policy}@*)
    PerspectiveProjection* PerspectiveProjection::clone() const
    {
        return new PerspectiveProjection(*this);
    }

    //(*@\Green{// Implementation of orthogonal projection matrices.}@*)

    //(*@\Green{// special constructor}@*)
    OrthogonalProjection::OrthogonalProjection(
            GLfloat aspect,
            GLfloat x_min, GLfloat x_max,
            GLfloat y_min, GLfloat y_max,
            GLfloat z_min, GLfloat z_max):
        _aspect(aspect),
        _x_min(x_min), _x_max(x_max),
        _y_min(y_min), _y_max(y_max),
        _z_min(z_min), _z_max(z_max)
    {
        _matrix[ 0] =  2.0f / (_x_max - _x_min) / _aspect;
        _matrix[ 5] =  2.0f / (_y_max - _y_min);
        _matrix[10] = -2.0f / (_z_max - _z_min);
        _matrix[12] = (_x_min + _x_max) / (_x_min - _x_max);
        _matrix[13] = (_y_min + _y_max) / (_y_min - _y_max);
        _matrix[14] = (_z_min + _z_max) / (_z_min - _z_max);
    }

    //(*@\Green{// setters}@*)
    GLvoid OrthogonalProjection::setAspectRatio(GLfloat aspect)
    {
        _aspect    = aspect;
        _matrix[0] = 2.0f / (_x_max - _x_min) / _aspect;
    }

    GLvoid OrthogonalProjection::setXMin(GLfloat x_min)
    {
        _x_min      = x_min;
        _matrix[ 0] = 2.0f / (_x_max - _x_min) / _aspect;
        _matrix[12] = (_x_min + _x_max) / (_x_min - _x_max);
    }

    GLvoid OrthogonalProjection::setXMax(GLfloat x_max)
    {
        _x_max      = x_max;
        _matrix[ 0] = 2.0f / (_x_max - _x_min) / _aspect;
        _matrix[12] = (_x_min + _x_max) / (_x_min - _x_max);
    }

    GLvoid OrthogonalProjection::setYMin(GLfloat y_min)
    {
        _y_min      = y_min;
        _matrix[ 5] = 2.0f / (_y_max - _y_min);
        _matrix[13] = (_y_min + _y_max) / (_y_min - _y_max);
    }

    GLvoid OrthogonalProjection::setYMax(GLfloat y_max)
    {
        _y_max      = y_max;
        _matrix[ 5] = 2.0f / (_y_max - _y_min);
        _matrix[13] = (_y_min + _y_max) / (_y_min - _y_max);
    }

    GLvoid OrthogonalProjection::setZMin(GLfloat z_min)
    {
        _z_min       = z_min;
        _matrix[10] = 2.0f / (_z_max - _z_min);
        _matrix[14] = (_z_min + _z_max) / (_z_min - _z_max);
    }

    GLvoid OrthogonalProjection::setZMax(GLfloat z_max)
    {
        _z_max      = z_max;
        _matrix[10] = 2.0f / (_z_max - _z_min);
        _matrix[14] = (_z_min + _z_max) / (_z_min - _z_max);
    }

    //(*@\Green{// redefined clone function required by smart pointers based on the deep copy ownership policy}@*)
    OrthogonalProjection* OrthogonalProjection::clone() const
    {
        return new OrthogonalProjection(*this);
    }

    //(*@\Green{// Implementation of world coordinate transformations.}@*)

    //(*@\Green{// special constructor}@*)
    LookAt::LookAt(const Cartesian3 &eye, const Cartesian3 &center, const Cartesian3 &up):
        _eye(eye), _center(center), _up(up)
    {
        Cartesian3 look(_eye);
        look -= _center;

        look.normalize();

        Cartesian3 right(_up);
        right ^= look;
        right.normalize();

        _up = look;
        _up ^= right;

        //(*@\Green{// first column}@*)
        _matrix[0] = (GLfloat)right[0];
        _matrix[1] = (GLfloat)  _up[0];
        _matrix[2] = (GLfloat) look[0];
        _matrix[3] = 0.0f;

        //(*@\Green{// second column}@*)
        _matrix[4] = (GLfloat)right[1];
        _matrix[5] = (GLfloat)  _up[1];
        _matrix[6] = (GLfloat) look[1];
        _matrix[7] = 0.0f;

        //(*@\Green{// third column}@*)
        _matrix[ 8] = (GLfloat)right[2];
        _matrix[ 9] = (GLfloat)  _up[2];
        _matrix[10] = (GLfloat) look[2];
        _matrix[11] = 0.0f;

        //(*@\Green{// fourth column}@*)
        _matrix[12] = -(GLfloat)(right * _eye);
        _matrix[13] = -(GLfloat)(  _up * _eye);
        _matrix[14] = -(GLfloat)( look * _eye);
        _matrix[15] = 1.0f;
    }

    //(*@\Green{// redefined clone function required by smart pointers based on the deep copy ownership policy}@*)
    LookAt* LookAt::clone() const
    {
        return new LookAt(*this);
    }
}
