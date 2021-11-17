//----------------------------------------------------------------------------------
// File:        Core/Math/SpecialGLTransformations.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef SPECIALGLTRANSFORMATIONS_H
#define SPECIALGLTRANSFORMATIONS_H

#include "GenericGLTransformations.h"
#include "Constants.h"

namespace cagd
{
    //(*@\Green{// Can be used to define translation matrices.}@*)
    class Translate: public GLTransformation
    {
    public:
        //(*@\Green{// default/special constructor}@*)
        Translate(GLfloat ux = 0.0f, GLfloat uy = 0.0f, GLfloat uz = 0.0f);

        //(*@\Green{// special constructor, the given direction vector will be normalized}@*)
        Translate(Cartesian3 &direction, GLfloat distance);

        //(*@\Green{// setters}@*)
        GLvoid setXDirectionalUnits(GLfloat ux);
        GLvoid setYDirectionalUnits(GLfloat uy);
        GLvoid setZDirectionalUnits(GLfloat uz);

        //(*@\Green{// redeclared clone function required by smart pointers based on the deep copy ownership policy}@*)
        Translate* clone() const;
    };

    //(*@\Green{// Can be used to define scaling transformations.}@*)
    class Scale: public GLTransformation
    {
    public:
        //(*@\Green{// default/special constructor}@*)
        Scale(GLfloat sx = 1.0f, GLfloat sy = 1.0f, GLfloat sz = 1.0f);

        //(*@\Green{// sets the scaling factors along axes $x$, $y$ and $z$}@*)
        GLvoid setScalingFactors(GLfloat sx, GLfloat sy, GLfloat sz);

        //(*@\Green{// redeclared clone function required by smart pointers based on the deep copy ownership policy}@*)
        Scale* clone() const;
    };

    //(*@\Green{// Can be used to define rotation matrices.}@*)
    class Rotate: public GLTransformation
    {
    protected:
        Cartesian3 _direction;
        GLfloat    _angle_in_radians;

    public:
        //(*@\Green{// default/special constructor}@*)
        Rotate(const Cartesian3 &direction = Cartesian3(1.0, 0.0, 0.0),
               GLfloat angle_in_radians = 0.0f);

        //(*@\Green{// sets the unit direction vector of the rotation axis}@*)
        GLvoid setDirection(const Cartesian3 &direction);

        //(*@\Green{// sets the rotation angle in radians}@*)
        GLvoid setAngle(GLfloat angle_in_radians);

        //(*@\Green{// redeclared clone function required by smart pointers based on the deep copy ownership policy}@*)
        Rotate* clone() const;
    };

    //(*@\Green{// Can be used to define perspective projection matrices.}@*)
    class PerspectiveProjection: public GLTransformation
    {
    protected:
        GLfloat _aspect;
        GLfloat _fov;
        GLfloat _f;
        GLfloat _near;
        GLfloat _far;

    public:
        //(*@\Green{// special constructor}@*)
        PerspectiveProjection(
                GLfloat aspect,
                GLfloat fov = 45.0f * DEG_TO_RADIAN,
                GLfloat near = 1.0f, GLfloat far = 1000.0f);

        //(*@\Green{// setters}@*)
        GLvoid setAspectRatio(GLfloat aspect);
        GLvoid setFieldOfView(GLfloat fov);
        GLvoid setNearClippingPlaneDistance(GLfloat near);
        GLvoid setFarClippingPlaneDistance(GLfloat far);

        //(*@\Green{// redeclared clone function required by smart pointers based on the deep copy ownership policy}@*)
        PerspectiveProjection* clone() const;
    };

    //(*@\Green{// Can be used to define orthogonal projection matrices.}@*)
    class OrthogonalProjection: public GLTransformation
    {
    protected:
        GLfloat _aspect;
        GLfloat _x_min, _x_max;
        GLfloat _y_min, _y_max;
        GLfloat _z_min, _z_max;

    public:
        //(*@\Green{// special constructor}@*)
        OrthogonalProjection(
            GLfloat aspect,
            GLfloat x_min, GLfloat x_max,
            GLfloat y_min, GLfloat y_max,
            GLfloat z_min, GLfloat z_max);

        //(*@\Green{// setters}@*)
        GLvoid setAspectRatio(GLfloat aspect);

        GLvoid setXMin(GLfloat x_min);
        GLvoid setXMax(GLfloat x_max);

        GLvoid setYMin(GLfloat y_min);
        GLvoid setYMax(GLfloat y_max);

        GLvoid setZMin(GLfloat z_min);
        GLvoid setZMax(GLfloat z_max);

        //(*@\Green{// redeclared clone function required by smart pointers based on the deep copy ownership policy}@*)
        OrthogonalProjection* clone() const;
    };

    //(*@\Green{// Can be used to define world coordinate transformations.}@*)
    class LookAt: public GLTransformation
    {
    protected:
        Cartesian3 _eye;
        Cartesian3 _center;
        Cartesian3 _up;

    public:
        //(*@\Green{// special constructor}@*)
        LookAt(const Cartesian3 &eye, const Cartesian3 &center, const Cartesian3 &up);

        //(*@\Green{// redeclared clone function required by smart pointers based on the deep copy ownership policy}@*)
        LookAt* clone() const;
    };
}

#endif //(*@\Green{// SPECIALGLTRANSFORMATIONS\_H}@*)
