//----------------------------------------------------------------------------------
// File:        Core/Geometry/Surfaces/Lights.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "Lights.h"
#include "../../Exceptions.h"

#include <new>

using namespace std;

namespace cagd
{
    //(*@\Green{// special constructor}@*)
    DirectionalLight::DirectionalLight(
        const Homogeneous3 &position,
        const Homogeneous3 &half_vector,
        const Color4       &ambient_intensity,
        const Color4       &diffuse_intensity,
        const Color4       &specular_intensity):

        _position(position),
        _half_vector(half_vector),
        _ambient_intensity(ambient_intensity),
        _diffuse_intensity(diffuse_intensity),
        _specular_intensity(specular_intensity)
    {
    }

    //(*@\Green{// setters}@*)
    GLvoid DirectionalLight::setPosition(const Homogeneous3 &position)
    {
        _position = position;
    }

    GLvoid DirectionalLight::setHalfVector(const Homogeneous3 &half_vector)
    {
        _half_vector = half_vector;
    }

    GLvoid DirectionalLight::setAmbientIntensity(const Color4 &ambient_intensity)
    {
        _ambient_intensity = ambient_intensity;
    }

    GLvoid DirectionalLight::setDiffuseIntensity(const Color4 &diffuse_intensity)
    {
        _diffuse_intensity = diffuse_intensity;
    }

    GLvoid DirectionalLight::setSpecularIntensity(const Color4 &specular_intensity)
    {
        _specular_intensity = specular_intensity;
    }

    //(*@\Green{// getters}@*)
    const GLfloat* DirectionalLight::addressOfPosition() const
    {
        return _position.address();
    }

    const GLfloat* DirectionalLight::addressOfHalfVector() const
    {
        return _half_vector.address();
    }

    const GLfloat* DirectionalLight::addressOfAmbientIntensity() const
    {
        return _ambient_intensity.address();
    }

    const GLfloat* DirectionalLight::addressOfDiffuseIntensity() const
    {
        return _diffuse_intensity.address();
    }

    const GLfloat* DirectionalLight::addressOfSpecularIntensity() const
    {
        return _specular_intensity.address();
    }

    //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    DirectionalLight* DirectionalLight::clone() const
    {
        return new (nothrow) DirectionalLight(*this);
    }

    //(*@\Green{// virtual default destructor}@*)
    DirectionalLight::~DirectionalLight() = default;

    //(*@\Green{// special constructor}@*)
    PointLight::PointLight(
        const Homogeneous3  &position,
        const Color4        &ambient_intensity,
        const Color4        &diffuse_intensity,
        const Color4        &specular_intensity,
        GLfloat             constant_attenuation,
        GLfloat             linear_attenuation,
        GLfloat             quadratic_attenuation):

        DirectionalLight(position, Homogeneous3(),
                         ambient_intensity, diffuse_intensity, specular_intensity),

        _constant_attenuation(constant_attenuation),
        _linear_attenuation(linear_attenuation),
        _quadratic_attenuation(quadratic_attenuation)
    {
        if (position.w() == 0.0f)
        {
            throw Exception("PointLight::PointLight - Wrong position.");
        }
    }

    //(*@\Green{// setters}@*)
    GLvoid PointLight::setConstantAttenuation(GLfloat constant_attenuation)
    {
        _constant_attenuation = constant_attenuation;
    }

    GLvoid PointLight::setLinearAttenuation(GLfloat linear_attenuation)
    {
        _linear_attenuation = linear_attenuation;
    }

    GLvoid PointLight::setQuadraticAttenuation(GLfloat quadratic_attenuation)
    {
        _quadratic_attenuation = quadratic_attenuation;
    }

    //(*@\Green{// getters}@*)
    GLfloat PointLight::constantAttenuation() const
    {
        return _constant_attenuation;
    }

    GLfloat PointLight::linearAttenuation() const
    {
        return _linear_attenuation;
    }

    GLfloat PointLight::quadraticAttenuation() const
    {
        return _quadratic_attenuation;
    }

    //(*@\Green{// redefined clone function required by smart pointers based on the deep copy ownership policy}@*)
    PointLight* PointLight::clone() const
    {
        return new (nothrow) PointLight(*this);
    }

    //(*@\Green{// virtual default destructor}@*)
    PointLight::~PointLight() = default;

    //(*@\Green{// special constructor}@*)
    Spotlight::Spotlight(
        const Homogeneous3 &position,
        const Color4       &ambient_intensity,
        const Color4       &diffuse_intensity,
        const Color4       &specular_intensity,
        GLfloat            constant_attenuation,
        GLfloat            linear_attenuation,
        GLfloat            quadratic_attenuation,
        const Cartesian3   &spot_direction,
        GLfloat            spot_cos_cutoff,
        GLfloat            spot_exponent):

        PointLight(position,
                   ambient_intensity, diffuse_intensity, specular_intensity,
                   constant_attenuation, linear_attenuation, quadratic_attenuation),

        _spot_direction(spot_direction),
        _spot_cos_cutoff(spot_cos_cutoff),
        _spot_exponent(spot_exponent)
    {
        if (position.w() == 0.0f)
        {
            throw Exception("Spotlight::Spotlight - Wrong position.");
        }

        if (_spot_cos_cutoff > 90.0f)
        {
            throw Exception("Spotlight::Spotlight - Wrong spot cosine cutoff.");
        }
    }

    //(*@\Green{// setters}@*)
    GLvoid Spotlight::setSpotDirection(const Cartesian3 &spot_direction)
    {
        _spot_direction = spot_direction;
    }

    GLvoid Spotlight::setSpotCosCutoff(GLfloat spot_cos_cutoff)
    {
        _spot_cos_cutoff = spot_cos_cutoff;
    }

    GLvoid Spotlight::setSpotExponent(GLfloat spot_exponent)
    {
        _spot_exponent = spot_exponent;
    }

    //(*@\Green{// getters}@*)
    const GLdouble* Spotlight::addressOfSpotDirection() const
    {
        return _spot_direction.address();
    }

    GLfloat Spotlight::spotCosCutoff() const
    {
        return _spot_cos_cutoff;
    }

    GLfloat Spotlight::spotExponent() const
    {
        return _spot_exponent;
    }

    //(*@\Green{// redefined clone function required by smart pointers based on the deep copy ownership policy}@*)
    Spotlight* Spotlight::clone() const
    {
        return new (nothrow) Spotlight(*this);
    }

    //(*@\Green{// virtual default destructor}@*)
    Spotlight::~Spotlight() = default;
}
