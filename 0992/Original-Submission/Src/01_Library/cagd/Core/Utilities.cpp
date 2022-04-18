//----------------------------------------------------------------------------------
// File:        Core/Utilities.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "Utilities.h"

#include <omp.h>
#include <GL/glew.h>

namespace cagd
{
    bool platformIsSupported()
    {
        return (omp_get_max_threads() > 1 && glewIsSupported("GL_VERSION_3_0"));
    }

    //(*@\Green{// {\color[rgb]{0, 0, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.0645161, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.129032, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.193548, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.258065, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.322581, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.387097, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.451613, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.516129, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.580645, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.645161, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.709677, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.774194, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.83871, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.903226, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 0.967742, 1}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.967742}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.903226}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.83871}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.774194}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.709677}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.645161}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.580645}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.516129}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.451613}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.387097}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.322581}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.258065}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.193548}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.129032}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0.0645162}\rule{0.2cm}{0.2cm}}{\color[rgb]{0, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.0645161, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.129032, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.193548, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.258065, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.322581, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.387097, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.451613, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.516129, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.580645, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.645161, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.709677, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.774194, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.83871, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.903226, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{0.967742, 1, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.967742, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.903226, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.83871, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.774194, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.709677, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.645161, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.580645, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.516129, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.451613, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.387097, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.322581, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.258065, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.193548, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.129032, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0.0645163, 0}\rule{0.2cm}{0.2cm}}{\color[rgb]{1, 0, 0}\rule{0.2cm}{0.2cm}}}\label{src:Utilities:coldToHotColormap:rainbow}@*)
    Color4 coldToHotColormap(GLfloat value, GLfloat min_value, GLfloat max_value) //(*@\label{src:Utilities:coldToHotColormap:start}@*)
    {
        Color4 color(1.0, 1.0, 1.0);

        if (value < min_value)
        {
           value = min_value;
        }

        if (value > max_value)
        {
           value = max_value;
        }

        float dv = max_value - min_value;

        if (value < (min_value + 0.25f * dv))
        {
           color.r() = 0.0;
           color.g() = 4.0f * (value - min_value) / dv;
        }
        else
        {
            if (value < (min_value + 0.5f * dv))
            {
               color.r() = 0.0f;
               color.b() = 1.0f + 4.0f * (min_value + 0.25f * dv - value) / dv;
            }
            else
            {
                if (value < (min_value + 0.75f * dv))
                {
                   color.r() = 4.0f * (value - min_value - 0.5f * dv) / dv;
                   color.b() = 0.0f;
                }
                else
                {
                   color.g() = 1.0f + 4.0f * (min_value + 0.75f * dv - value) / dv;
                   color.b() = 0.0f;
                }
            }
        }

        return color;
    } //(*@\label{src:Utilities:coldToHotColormap:end}@*)

    const char* openGLTypeToString(GLenum type)
    {
        switch (type)
        {
        case GL_FLOAT:                                     return "float";
        case GL_FLOAT_VEC2:                                return "vec2";
        case GL_FLOAT_VEC3:                                return "vec3";
        case GL_FLOAT_VEC4:                                return "vec4";
        case GL_INT:                                       return "int";
        case GL_INT_VEC2:                                  return "ivec2";
        case GL_INT_VEC3:                                  return "ivec3";
        case GL_INT_VEC4:                                  return "ivec4";
        case GL_UNSIGNED_INT:                              return "unsigned int";
        case GL_UNSIGNED_INT_VEC2:                         return "uvec2";
        case GL_UNSIGNED_INT_VEC3:                         return "uvec3";
        case GL_UNSIGNED_INT_VEC4:                         return "uvec4";
        case GL_BOOL:                                      return "bool";
        case GL_BOOL_VEC2:                                 return "bvec2";
        case GL_BOOL_VEC3:                                 return "bvec3";
        case GL_BOOL_VEC4:                                 return "bvec4";
        case GL_FLOAT_MAT2:                                return "mat2";
        case GL_FLOAT_MAT3:                                return "mat3";
        case GL_FLOAT_MAT4:                                return "mat4";
        case GL_FLOAT_MAT2x3:                              return "mat2x3";
        case GL_FLOAT_MAT2x4:                              return "mat2x4";
        case GL_FLOAT_MAT3x2:                              return "mat3x2";
        case GL_FLOAT_MAT3x4:                              return "mat3x4";
        case GL_FLOAT_MAT4x2:                              return "mat4x2";
        case GL_FLOAT_MAT4x3:                              return "mat4x3";
        case GL_SAMPLER_1D:                                return "sampler1D";
        case GL_SAMPLER_2D:                                return "sampler2D";
        case GL_SAMPLER_3D:                                return "sampler3D";
        case GL_SAMPLER_CUBE:                              return "samplerCube";
        case GL_SAMPLER_1D_SHADOW:                         return "sampler1DShadow";
        case GL_SAMPLER_2D_SHADOW:                         return "sampler2DShadow";
        case GL_SAMPLER_1D_ARRAY:                          return "sampler1DArray";
        case GL_SAMPLER_2D_ARRAY:                          return "sampler2DArray";
        case GL_SAMPLER_1D_ARRAY_SHADOW:                   return "sampler1DArrayShadow";
        case GL_SAMPLER_2D_ARRAY_SHADOW:                   return "sampler2DArrayShadow";
        case GL_SAMPLER_2D_MULTISAMPLE:                    return "sampler2DMS";
        case GL_SAMPLER_2D_MULTISAMPLE_ARRAY:              return "sampler2DMSArray";
        case GL_SAMPLER_CUBE_SHADOW:                       return "samplerCubeShadow";
        case GL_SAMPLER_BUFFER:                            return "samplerBuffer";
        case GL_SAMPLER_2D_RECT:                           return "sampler2DRect";
        case GL_SAMPLER_2D_RECT_SHADOW:                    return "sampler2DRectShadow";
        case GL_INT_SAMPLER_1D:                            return "isampler1D";
        case GL_INT_SAMPLER_2D:                            return "isampler2D";
        case GL_INT_SAMPLER_3D:                            return "isampler3D";
        case GL_INT_SAMPLER_CUBE:                          return "isamplerCube";
        case GL_INT_SAMPLER_1D_ARRAY:                      return "isampler1DArray";
        case GL_INT_SAMPLER_2D_ARRAY:                      return "isampler2DArray";
        case GL_INT_SAMPLER_2D_MULTISAMPLE:                return "isampler2DMS";
        case GL_INT_SAMPLER_2D_MULTISAMPLE_ARRAY:          return "isampler2DMSArray";
        case GL_INT_SAMPLER_BUFFER:                        return "isamplerBuffer";
        case GL_INT_SAMPLER_2D_RECT:                       return "isampler2DRect";
        case GL_UNSIGNED_INT_SAMPLER_1D:                   return "usampler1D";
        case GL_UNSIGNED_INT_SAMPLER_2D:                   return "usampler2D";
        case GL_UNSIGNED_INT_SAMPLER_3D:                   return "usampler3D";
        case GL_UNSIGNED_INT_SAMPLER_CUBE:                 return "usamplerCube";
        case GL_UNSIGNED_INT_SAMPLER_1D_ARRAY:             return "usampler2DArray";
        case GL_UNSIGNED_INT_SAMPLER_2D_ARRAY:             return "usampler2DArray";
        case GL_UNSIGNED_INT_SAMPLER_2D_MULTISAMPLE:       return "usampler2DMS";
        case GL_UNSIGNED_INT_SAMPLER_2D_MULTISAMPLE_ARRAY: return "usampler2DMSArray";
        case GL_UNSIGNED_INT_SAMPLER_BUFFER:               return "usamplerBuffer";
        case GL_UNSIGNED_INT_SAMPLER_2D_RECT:              return "usampler2DRect";
        default:                                           return "other";
        }
    }
}
