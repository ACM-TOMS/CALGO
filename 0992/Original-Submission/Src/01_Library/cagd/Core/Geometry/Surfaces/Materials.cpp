//----------------------------------------------------------------------------------
// File:        Core/Geometry/Surfaces/Materials.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "Materials.h"

#include <new>

using namespace std;

namespace cagd
{
    //(*@\Green{// default/special constructor}@*)
    Material::Material(const Color4 &ambient,  const Color4 &diffuse,
                       const Color4 &specular, const Color4 &emission,
                       GLfloat shininess):
        _ambient(ambient),
        _diffuse(diffuse),
        _specular(specular),
        _emission(emission),
        _shininess(shininess)
    {
    }

    //(*@\Green{// setters}@*)
    GLvoid Material::setAmbientReflectionCoefficients(
                GLfloat r, GLfloat g, GLfloat b, GLfloat a)
    {
        _ambient.r() = r;
        _ambient.g() = g;
        _ambient.b() = b;
        _ambient.a() = a;
    }

    GLvoid Material::setDiffuseReflectionCoefficients(
                GLfloat r, GLfloat g, GLfloat b, GLfloat a)
    {
        _diffuse.r() = r;
        _diffuse.g() = g;
        _diffuse.b() = b;
        _diffuse.a() = a;
    }

    GLvoid Material::setSpecularReflectionCoefficients(
                GLfloat r, GLfloat g, GLfloat b, GLfloat a)
    {
        _specular.r() = r;
        _specular.g() = g;
        _specular.b() = b;
        _specular.a() = a;
    }

    GLvoid Material::setEmissionColor(
                GLfloat r, GLfloat g, GLfloat b, GLfloat a)
    {
        _emission.r() = r;
        _emission.g() = g;
        _emission.b() = b;
        _emission.a() = a;
    }

    GLvoid Material::setAmbientReflectionCoefficients(const Color4 &c)
    {
        _ambient = c;
    }

    GLvoid Material::setDiffuseReflectionCoefficients(const Color4 &c)
    {
        _diffuse = c;
    }

    GLvoid Material::setSpecularReflectionCoefficients(const Color4 &c)
    {
        _specular = c;
    }

    GLvoid Material::setEmissionColor(const Color4 &c)
    {
        _emission = c;
    }

    GLvoid Material::setShininess(GLfloat shininess)
    {
        _shininess = shininess;
    }

    GLvoid Material::setTransparency(GLfloat alpha)
    {
        _ambient.a() = _diffuse.a() = _specular.a() = alpha;
    }

    //(*@\Green{// getters}@*)
    const GLfloat* Material::addressOfAmbientReflectionCoefficients() const
    {
        return _ambient.address();
    }

    const GLfloat* Material::addressOfDiffuseReflectionCoefficients() const
    {
        return _diffuse.address();
    }

    const GLfloat* Material::addressOfSpecularReflectionCoefficients() const
    {
        return _specular.address();
    }

    const GLfloat* Material::addressOfEmissionColor() const
    {
        return _emission.address();
    }

    GLfloat Material::shininess() const
    {
        return _shininess;
    }

    GLboolean Material::isTransparent() const
    {
        return (_diffuse.a() < 1.0f);
    }

    //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    Material* Material::clone() const
    {
        return new (nothrow) Material(*this);
    }

    //(*@\Green{// predefined materials}@*)
    namespace materials //(*@\label{src:predefined_materials:start}@*)
    {
        Material black_plastique = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/black_plastique.pdf}}@*)
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.4f),
                        Color4(0.020000f, 0.020000f, 0.020000f, 0.6f),
                        Color4(0.600000f, 0.600000f, 0.600000f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        32.0f),

                 black_rubber = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/black_rubber.pdf}}@*)
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.4f),
                        Color4(0.010000f, 0.010000f, 0.010000f, 0.6f),
                        Color4(0.500000f, 0.500000f, 0.500000f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        32.0f),

                 brass = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/brass.pdf}}@*)
                        Color4(0.329412f, 0.223529f, 0.027451f, 0.4f),
                        Color4(0.780392f, 0.568627f, 0.113725f, 0.6f),
                        Color4(0.992157f, 0.941176f, 0.807843f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        27.8974f),

                 bronze = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/bronze.pdf}}@*)
                        Color4(0.250000f, 0.148000f, 0.064750f, 0.4f),
                        Color4(0.400000f, 0.236800f, 0.103600f, 0.6f),
                        Color4(0.774597f, 0.458561f, 0.200621f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        76.8f),

                 chrome = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/chrome.pdf}}@*)
                        Color4(0.250000f, 0.250000f, 0.250000f, 0.4f),
                        Color4(0.400000f, 0.400000f, 0.400000f, 0.6f),
                        Color4(0.774597f, 0.774597f, 0.774597f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        76.8f),

                 copper = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/copper.pdf}}@*)
                        Color4(0.191250f, 0.073500f, 0.022500f, 0.4f),
                        Color4(0.703800f, 0.270480f, 0.082800f, 0.6f),
                        Color4(0.256777f, 0.137622f, 0.086014f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        51.2f),

                 emerald = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/emerald.pdf}}@*)
                        Color4(0.021500f, 0.174500f, 0.021500f, 0.4f),
                        Color4(0.075680f, 0.614240f, 0.075680f, 0.6f),
                        Color4(0.633000f, 0.727811f, 0.633000f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        76.8f),

                 gold = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/gold.pdf}}@*)
                        Color4(0.247250f, 0.224500f, 0.064500f, 0.4f),
                        Color4(0.346150f, 0.314300f, 0.090300f, 0.6f),
                        Color4(0.797357f, 0.723991f, 0.208006f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        83.2f),

                 jade = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/jade.pdf}}@*)
                        Color4(0.135000f, 0.222500f, 0.157500f, 0.4f),
                        Color4(0.540000f, 0.890000f, 0.630000f, 0.6f),
                        Color4(0.316228f, 0.316228f, 0.316228f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        12.8f),

                 obsidian = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/obsidian.pdf}}@*)
                        Color4(0.053750f, 0.050000f, 0.066250f, 0.4f),
                        Color4(0.182750f, 0.170000f, 0.225250f, 0.6f),
                        Color4(0.332741f, 0.328634f, 0.346435f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        38.4f),

                 pearl = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/pearl.pdf}}@*)
                        Color4(0.250000f, 0.207250f, 0.207250f, 0.4f),
                        Color4(1.000000f, 0.829000f, 0.829000f, 0.6f),
                        Color4(0.296648f, 0.296648f, 0.296648f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        11.264f),

                 pewter = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/pewter.pdf}}@*)
                        Color4(0.105882f, 0.058824f, 0.113725f, 0.4f),
                        Color4(0.427451f, 0.470588f, 0.541176f, 0.6f),
                        Color4(0.333333f, 0.333333f, 0.521569f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        9.84615f),

                 polished_bronze = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/polished_bronze.pdf}}@*)
                        Color4(0.212500f, 0.127500f, 0.054000f, 0.4f),
                        Color4(0.714000f, 0.428400f, 0.181440f, 0.6f),
                        Color4(0.393548f, 0.271906f, 0.166721f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        25.6f),

                 polished_copper = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/polished_copper.pdf}}@*)
                        Color4(0.229500f, 0.088250f, 0.0275000f, 0.4f),
                        Color4(0.550800f, 0.211800f, 0.0660000f, 0.6f),
                        Color4(0.580594f, 0.223257f, 0.0695701f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.0000000f, 0.0f),
                        12.8f),

                 polished_gold = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/polished_gold.pdf}}@*)
                        Color4(0.247250f, 0.199500f, 0.074500f, 0.4f),
                        Color4(0.751640f, 0.606480f, 0.226480f, 0.6f),
                        Color4(0.628281f, 0.555802f, 0.366065f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        51.2f),

                 polished_silver = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/polished_silver.pdf}}@*)
                        Color4(0.192250f, 0.192250f, 0.192250f, 0.4f),
                        Color4(0.507540f, 0.507540f, 0.507540f, 0.6f),
                        Color4(0.508273f, 0.508273f, 0.508273f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        51.2f),

                 ruby = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/ruby.pdf}}@*)
                        Color4(0.174500f, 0.011750f, 0.011750f, 0.4f),
                        Color4(0.614240f, 0.041360f, 0.041360f, 0.6f),
                        Color4(0.727811f, 0.626959f, 0.626959f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        76.8f),

                 silver = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{\includegraphics{pdf/silver.pdf}}@*)
                        Color4(0.231250f, 0.231250f, 0.231250f, 0.4f),
                        Color4(0.277500f, 0.277500f, 0.277500f, 0.6f),
                        Color4(0.773911f, 0.773911f, 0.773911f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        89.6f),

                turquoise = Material(//(*@\hspace*{\fill}\makebox[4.65cm][l]{$\includegraphics{pdf/turquoise.pdf}^{\footnote[9]{\Green{The model file mouse.off is the creation of the author. It was rendered by means of the proposed function library.}}}$}@*)
                        Color4(0.100000f, 0.187250f, 0.174500f, 0.4f),
                        Color4(0.396000f, 0.741510f, 0.691020f, 0.6f),
                        Color4(0.297254f, 0.308290f, 0.306678f, 0.8f),
                        Color4(0.000000f, 0.000000f, 0.000000f, 0.0f),
                        12.8f);
    }//(*@\label{src:predefined_materials:end}@*)
}
