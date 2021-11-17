//----------------------------------------------------------------------------------
// File:        Core/Geometry/Surfaces/Materials.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef MATERIALS_H
#define MATERIALS_H

#include "../Coordinates/Colors4.h"

namespace cagd
{
    class Material
    {
    private:
        //(*@\Green{// ambient, diffuse, specular reflection coefficients and emissive color components}@*)
        Color4	_ambient, _diffuse, _specular, _emission;

        //(*@\Green{// shininess of the material, its value should be in the interval $\left[0,128\right]$}@*)
        GLfloat	_shininess;

    public:
        //(*@\Green{// default/special constructor}@*)
        Material(const Color4 &ambient	= Color4(),
                 const Color4 &diffuse	= Color4(),
                 const Color4 &specular = Color4(),
                 const Color4 &emission = Color4(),
                 GLfloat shininess = 128.0f);

        //(*@\Green{// setters}@*)
        GLvoid setAmbientReflectionCoefficients(const Color4 &c);
        GLvoid setAmbientReflectionCoefficients(
                    GLfloat r, GLfloat g, GLfloat b, GLfloat a = 1.0f);

        GLvoid setDiffuseReflectionCoefficients(const Color4 &c);
        GLvoid setDiffuseReflectionCoefficients(
                    GLfloat r, GLfloat g, GLfloat b, GLfloat a = 1.0f);

        GLvoid setSpecularReflectionCoefficients(const Color4 &c);
        GLvoid setSpecularReflectionCoefficients(
                    GLfloat r, GLfloat g, GLfloat b, GLfloat a = 1.0f);

        GLvoid setEmissionColor(const Color4 &c);
        GLvoid setEmissionColor(
                    GLfloat r, GLfloat g, GLfloat b, GLfloat a = 1.0f);

        GLvoid setShininess(GLfloat shininess);
        GLvoid setTransparency(GLfloat alpha);

        //(*@\Green{// getters}@*)
        const GLfloat* addressOfAmbientReflectionCoefficients()  const;
        const GLfloat* addressOfDiffuseReflectionCoefficients()  const;
        const GLfloat* addressOfSpecularReflectionCoefficients() const;
        const GLfloat* addressOfEmissionColor() const;

        GLfloat        shininess() const;
        GLboolean      isTransparent() const;

        //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        Material* clone() const;
    };

    //(*@\Green{// predefined materials}@*)
    namespace materials
    {
        extern Material black_plastique,
                        black_rubber,
                        brass,
                        bronze,
                        chrome,
                        copper,
                        emerald,
                        gold,
                        jade,
                        obsidian,
                        pearl,
                        pewter,
                        polished_bronze,
                        polished_copper,
                        polished_gold,
                        polished_silver,
                        ruby,
                        silver,
                        turquoise;
    }
}

#endif //(*@\Green{// MATERIALS\_H}@*)
