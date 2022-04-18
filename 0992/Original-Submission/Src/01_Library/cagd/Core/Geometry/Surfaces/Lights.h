//----------------------------------------------------------------------------------
// File:        Core/Geometry/Surfaces/Lights.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef LIGHTS_H
#define LIGHTS_H

#include "../Coordinates/Cartesians3.h"
#include "../Coordinates/Colors4.h"
#include "../Coordinates/Homogeneous3.h"

namespace cagd
{
    class DirectionalLight
    {
    protected:
        Homogeneous3 _position, _half_vector;
        Color4       _ambient_intensity, _diffuse_intensity, _specular_intensity;

    public:
        //(*@\Green{// special constructor}@*)
        DirectionalLight(
            const Homogeneous3 &position,
            const Homogeneous3 &half_vector,
            const Color4       &ambient_intensity,
            const Color4       &diffuse_intensity,
            const Color4       &specular_intensity);

        //(*@\Green{// setters}@*)
        GLvoid setPosition(const Homogeneous3 &position);
        GLvoid setHalfVector(const Homogeneous3 &half_vector);
        GLvoid setAmbientIntensity(const Color4 &ambient_intensity);
        GLvoid setDiffuseIntensity(const Color4 &diffuse_intensity);
        GLvoid setSpecularIntensity(const Color4 &specular_intensity);

        //(*@\Green{// getters}@*)
        const GLfloat* addressOfPosition() const;
        const GLfloat* addressOfHalfVector() const;
        const GLfloat* addressOfAmbientIntensity() const;
        const GLfloat* addressOfDiffuseIntensity() const;
        const GLfloat* addressOfSpecularIntensity() const;

        //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        virtual DirectionalLight* clone() const;

        //(*@\Green{// virtual destructor}@*)
        virtual ~DirectionalLight();
    };

    class PointLight: public DirectionalLight
    {
    protected:
        GLfloat _constant_attenuation,
                _linear_attenuation,
                _quadratic_attenuation;

    public:
        //(*@\Green{// special constructor}@*)
        PointLight(
            const Homogeneous3  &position,
            const Color4        &ambient_intensity,
            const Color4        &diffuse_intensity,
            const Color4        &specular_intensity,
            GLfloat             constant_attenuation,
            GLfloat             linear_attenuation,
            GLfloat             quadratic_attenuation);

        //(*@\Green{// setters}@*)
        GLvoid setConstantAttenuation(GLfloat constant_attenuation);
        GLvoid setLinearAttenuation(GLfloat linear_attenuation);
        GLvoid setQuadraticAttenuation(GLfloat quadratic_attenuation);

        //(*@\Green{// getters}@*)
        GLfloat constantAttenuation() const;
        GLfloat linearAttenuation() const;
        GLfloat quadraticAttenuation() const;

        //(*@\Green{// redeclared clone function required by smart pointers based on the deep copy ownership policy}@*)
        PointLight* clone() const;

        //(*@\Green{// virtual destructor}@*)
        virtual ~PointLight();
    };

    class Spotlight: public PointLight
    {
    private:
        Cartesian3 _spot_direction;
        GLfloat    _spot_cos_cutoff, _spot_exponent;

    public:
        //(*@\Green{// special constructor}@*)
        Spotlight(
            const Homogeneous3& position,
            const Color4        &ambient_intensity,
            const Color4        &diffuse_intensity,
            const Color4        &specular_intensity,
            GLfloat             constant_attenuation,
            GLfloat             linear_attenuation,
            GLfloat             quadratic_attenuation,
            const Cartesian3    &spot_direction,
            GLfloat             spot_cos_cutoff,
            GLfloat             spot_exponent);

        //(*@\Green{// setters}@*)
        GLvoid setSpotDirection(const Cartesian3 &spot_direction);
        GLvoid setSpotCosCutoff(GLfloat spot_cos_cutoff);
        GLvoid setSpotExponent(GLfloat spot_exponent);

        //(*@\Green{// getters}@*)
        const GLdouble* addressOfSpotDirection() const;
        GLfloat spotCosCutoff() const;
        GLfloat spotExponent() const;

        //(*@\Green{// redeclared clone function required by smart pointers based on the deep copy ownership policy}@*)
        Spotlight* clone() const;

        //(*@\Green{// virtual default destructor}@*)
        virtual ~Spotlight();
    };
}

#endif //(*@\Green{// LIGHTS\_H}@*)
