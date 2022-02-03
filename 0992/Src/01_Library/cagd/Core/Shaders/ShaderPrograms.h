//----------------------------------------------------------------------------------
// File:        Core/Shaders/ShaderPrograms.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef SHADERPROGRAMS_H
#define SHADERPROGRAMS_H

#include <GL/glew.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "../Geometry/Surfaces/Lights.h"
#include "../Geometry/Surfaces/Materials.h"

namespace cagd
{
    class ShaderProgram
    {
        //(*@\Green{// friend classes that will be discussed in Listings \mref{src:GenericCurves3.h}--\mref{src:GenericCurves3.cpp}, \mref{src:LinearCombinations3.h}--\mref{src:LinearCombinations3.cpp}, \mref{src:TriangleMeshes3.h}--\mref{src:TriangleMeshes3.cpp}}@*)
        //(*@\Green{// and \mref{src:TensorProductSurfaces3.h}--\mref{src:TensorProductSurfaces3.cpp}, respectively}@*)
        friend class GenericCurve3;
        friend class LinearCombination3;
        friend class TriangleMesh3;
        friend class TensorProductSurface3;

    public:
        //(*@\Green{// nested class that represent different types of shaders}@*)
        class Shader
        {
            friend class ShaderProgram;

        public:
            //(*@\Green{// possible shader types}@*)
            enum Type{VERTEX = 0,
                      FRAGMENT,
                      COMPUTE,
                      GEOMETRY,
                      TESS_CONTROL,
                      TESS_EVALUATION};

        private:
            GLuint      _id;
            GLint       _compile_status;
            Type        _type;
            std::string _source;

        public:
            //(*@\Green{// default constructor}@*)
            Shader();

            //(*@\Green{// copy constructor}@*)
            Shader(const Shader &shader);

            //(*@\Green{// assignment operator}@*)
            Shader& operator =(const Shader &rhs);

            //(*@\Green{// creates, loads and compiles from source file}@*)
            GLboolean createLoadAndCompileFromSourceFile(
                    Type type, const std::string &file_name,
                    GLboolean logging_is_enabled = GL_FALSE,
                    std::ostream &output = std::cout);

            //(*@\Green{// either deletes, or flags for deletion the shader}@*)
            GLvoid flagForDeletion() const;

            //(*@\Green{// destructor}@*)
            ~Shader();
        };

    public:
        //(*@\Green{// nested class that represents active uniform and attribute variables}@*)
        class ActiveVariable
        {
        private:
            GLenum  _type;
            GLint   _array_size;
            GLint   _location;

        public:
            //(*@\Green{// default/special constructor}@*)
            ActiveVariable(GLenum type = GL_NONE, GLint arraySize = 0,
                           GLint location = -1);

            //(*@\Green{// getters}@*)
            GLenum        type() const;
            GLint         arraySize() const;
            GLint         location()  const;
            const GLchar* typeDescription() const;
        };

        typedef ActiveVariable Attribute;
        typedef ActiveVariable Uniform;

    private:
        GLuint                           _id;
        GLint                            _link_status;
        std::vector<Shader>              _shader;
        std::map<std::string, Attribute> _attribute;
        std::map<std::string, Uniform>   _uniform;

        std::string                      _position_attribute_name;
        std::string                      _normal_attribute_name;
        std::string                      _color_attribute_name;
        std::string                      _texture_attribute_name;

    public:
        //(*@\Green{// default constructor}@*)
        //(*@\Green{//}@*)
        //(*@\Green{// Creates an empty unlinked shader program and sets the names of position, normal,}@*)
        //(*@\Green{// color and texture attributes to "position", "normal", "color" and "texture",}@*)
        //(*@\Green{// respectively.}@*)
        //(*@\Green{//}@*)
        //(*@\Green{// If required, these attribute names can be changed by using one of the methods}@*)
        //(*@\Green{// setPositionAttributeName, setNormalAttributeName, setColorAttributeName or}@*)
        //(*@\Green{// setTextureAttributeName.}@*)
        //(*@\Green{//}@*)
        //(*@\Green{// By default, it is assumed that, in the vertex shader, these attributes are defined as}@*)
        //(*@\Green{// the input variables:}@*)
        //(*@\Green{//}@*)
        //(*@\Green{// in vec3 position;}@*)
        //(*@\Green{// in vec3 normal;}@*)
        //(*@\Green{// in vec4 color;}@*)
        //(*@\Green{// in vec4 texture;}@*)
        //(*@\Green{//}@*)
        //(*@\Green{// The reason for this is that our rendering methods are based on vertex buffer objects}@*)
        //(*@\Green{// that describe:}@*)
        //(*@\Green{// \ \ - positions and unit normals as packages of 3 floats;}@*)
        //(*@\Green{// \ \ - (potentially projective) texture elements and colors as packages of 4 floats.}@*)
        //(*@\Green{//}@*)
        //(*@\Green{// The rendering methods GenericCurve3::renderDerivatives, LinearCombination3::renderData}@*)
        //(*@\Green{// and TensorProductSurface3::renderData require only the attribute "position", while the}@*)
        //(*@\Green{// constant color of the rendered derivatives, control polygons and control nets can be}@*)
        //(*@\Green{// specified as a uniform variable as it is illustrated in the next simple vertex and}@*)
        //(*@\Green{// fragment shaders:}@*)
        //(*@\Green{//}@*)
        //    ----------
        //    color.vert
        //    ----------
        //    uniform mat4 PVM;      //(*@\Green{// product of projection, view and model matrices}@*)
        //    in      vec3 position; //(*@\Green{// input attribute}@*)
        //
        //    void main()
        //    {
        //        gl_Position = PVM * vec4(position, 1.0);
        //    }
        //
        //    ----------
        //    color.frag
        //    ----------
        //    uniform vec4 color;           //(*@\Green{// constant/flat color}@*)
        //    out     vec4 fragment_color;  //(*@\Green{// output of the fragment shader}@*)
        //
        //    void main()
        //    {
        //        fragment_color = color;
        //    }
        //
        //(*@\Green{// However, the rendering method TriangleMesh3::render assumes that a successfully}@*)
        //(*@\Green{// compiled and linked shader program provides active attributes for positions,}@*)
        //(*@\Green{// unit normals and texture elements, i.e., in this case, one should use vertex and}@*)
        //(*@\Green{// fragment shaders similar to the next examples:}@*)
        //(*@\Green{//}@*)
        //    -----------
        //    shader.vert
        //    -----------
        //    uniform mat4 VM;  //(*@\Green{// product of view and model matrices}@*)
        //    uniform mat4 PVM; //(*@\Green{// product of projection, view and model matrices}@*)
        //    uniform mat4 N;   //(*@\Green{// normal matrix, i.e., transposed inverse of VM}@*)
        //
        //    //(*@\Green{// input attributes}@*)
        //    in      vec3 position;
        //    in      vec3 normal;
        //    in      vec4 color;
        //    in      vec4 texture;
        //
        //    //(*@\Green{// outputs of the vertex shader}@*)
        //    out     vec3 interpolated_position;
        //    out     vec3 interpolated_normal;
        //    out     vec4 interpolated_color;
        //    out     vec4 interpolated_texture;
        //
        //    void main()
        //    {
        //        //(*@\Green{// transform the normal vector to the eye space, then normalize it}@*)
        //        interpolated_normal = normalize(vec3(N * vec4(normal, 0.0)));
        //
        //        //(*@\Green{// transform the vertex position to the eye space}@*)
        //        //(*@\Green{// (it may be useful in case of point and spot lights)}@*)
        //        interpolated_position = vec3(VM * vec4(position, 1.0));
        //
        //        //(*@\Green{// interpolate texture coordinates and color components}@*)
        //        interpolated_color   = color;
        //        interpolated_texture = texture;
        //
        //        //(*@\Green{// convert the given vertex to clip coordinates and pass along}@*)
        //        gl_Position = PVM * vec4(position, 1.0);
        //    }
        //
        //    -----------
        //    shader.frag
        //    -----------
        //    in  vec3 interpolated_position;
        //    in  vec3 interpolated_normal;
        //    in  vec4 interpolated_color;
        //    in  vec4 interpolated_texture;
        //
        //    ...
        //
        //    out vec4 fragment_color; //(*@\Green{// output of the fragment shader}@*)
        //
        //    void main()
        //    {
        //        //(*@\Green{// one should use all interpolated values computed by the vertex shader,}@*)
        //        //(*@\Green{// otherwise the compiler will optimize out any inactive attributes and}@*)
        //        //(*@\Green{// the rendering method TriangleMesh3::render will return GL\_FALSE,}@*)
        //        //(*@\Green{// since it will be not able to find the locations of all required}@*)
        //        //(*@\Green{// active attributes that have to be passed to function calls like}@*)
        //        //(*@\Green{// glVertexAttribPointer}@*)
        //
        //        ...
        //
        //        fragment_color = ...;
        //    }
        ShaderProgram();

        //(*@\Green{// copy constructor}@*)
        ShaderProgram(const ShaderProgram &program);

        //(*@\Green{// assignment operator}@*)
        ShaderProgram& operator =(const ShaderProgram &rhs);

        //(*@\Green{// Tries to create, load and compile a shader from the specified source file.}@*)
        //(*@\Green{// In case of success, a new Shader object will be inserted at the end of the}@*)
        //(*@\Green{// std::vector$<$Shader$>$ \_shader.}@*)
        GLboolean attachNewShaderFromSourceFile(
                Shader::Type type, const std::string &file_name,
                GLboolean logging_is_enabled = GL_FALSE,
                std::ostream &output = std::cout);

        //(*@\Green{// Attaches an existing pre-compiled shader object to the shader program. The index}@*)
        //(*@\Green{// corresponds to the insertion order of the underlying shader. If the provided}@*)
        //(*@\Green{// index is greater than or equal to the size of the std::vector$<$Shader$>$ \_shader, the}@*)
        //(*@\Green{// method returns GL\_FALSE. If the specified shader exists and it is already attached}@*)
        //(*@\Green{// to the shader program, the method will also return GL\_FALSE.}@*)
        GLboolean attachExistingShader(GLuint index);

        //(*@\Green{// Detaches an existing shader object from the shader program. The index corresponds}@*)
        //(*@\Green{// to the insertion order of the underlying shader. In case of success, the shader}@*)
        //(*@\Green{// will be not deleted, it will simply be not considered in a later linking process}@*)
        //(*@\Green{// and, if required, it can be reattached later on. If the provided index is greater than}@*)
        //(*@\Green{// or equal to the size of the std::vector$<$Shader$>$ \_shader, the method returns GL\_FALSE.}@*)
        GLboolean detachExistingShader(GLuint index);

        //(*@\Green{// Links the attached pre-compiled shader objects. In case of success, the method}@*)
        //(*@\Green{// also loads the names, types, array sizes and locations of active uniforms and}@*)
        //(*@\Green{// attributes. The names and associated properties of these variables will}@*)
        //(*@\Green{// be stored in the private maps std::map$<$std::string, Uniform$>$ \_uniform and}@*)
        //(*@\Green{// std::map$<$std::string, Attribute$>$ \_attribute.}@*)
        GLboolean linkAttachedShaders(GLboolean logging_is_enabled = GL_FALSE,
                                      std::ostream &output = std::cout);

        //(*@\Green{// Assuming that the shader program was already compiled and successfully linked,}@*)
        //(*@\Green{// the method enables and validates its usage.}@*)
        GLboolean enable(GLboolean logging_is_enabled = GL_FALSE,
                         std::ostream &output = std::cout) const;

        //(*@\Green{// The method disables the underlying shader program.}@*)
        GLvoid    disable() const;

        //(*@\Green{// getters}@*)
        GLint     positionAttributeLocation() const;
        GLint     normalAttributeLocation() const;
        GLint     colorAttributeLocation() const;
        GLint     textureAttributeLocation() const;

        //(*@\Green{// Returns the properties of an active attribute that corresponds to the given name.}@*)
        //(*@\Green{// If there is no such active attribute, the default attribute object}@*)
        //(*@\Green{// will be returned -- the type, array size and location of which is set}@*)
        //(*@\Green{// to GL\_NONE/"other", 0 and -1, respectively.}@*)
        Attribute activeAttribute(const std::string &name) const;

        //(*@\Green{// Returns the properties of an active uniform that corresponds to the given name.}@*)
        //(*@\Green{// If there is no such active uniform, the default uniform object}@*)
        //(*@\Green{// will be returned -- the type, array size and location of which is set}@*)
        //(*@\Green{// to GL\_NONE/"other", 0 and -1, respectively.}@*)
        Uniform   activeUniform(const std::string &name) const;

        //(*@\Green{// setters}@*)
        GLvoid    setPositionAttributeName(const std::string &name); //(*@\label{src:ShaderProgram:setPositionAttributeName:declaration}@*)
        GLvoid    setNormalAttributeName(const std::string &name);
        GLvoid    setColorAttributeName(const std::string &name);
        GLvoid    setTextureAttributeName(const std::string &name);

        GLboolean setUniformValue1i(const std::string &name, GLint value) const; //(*@\label{src:ShaderProgram:setUniformVariables:declaration:start}@*)
        GLboolean setUniformValue2i(const std::string &name,
                                    GLint value_0, GLint value_1) const;
        GLboolean setUniformValue3i(const std::string &name,
                                    GLint value_0, GLint value_1, GLint value_2) const;
        GLboolean setUniformValue4i(const std::string &name,
                                    GLint value_0, GLint value_1,
                                    GLint value_2, GLint value_3) const;

        GLboolean setUniformValue1ui(const std::string &name, GLuint value) const;
        GLboolean setUniformValue2ui(const std::string &name,
                                     GLuint value_0, GLuint value_1) const;
        GLboolean setUniformValue3ui(const std::string &name,
                                     GLuint value_0, GLuint value_1,
                                     GLuint value_2) const;
        GLboolean setUniformValue4ui(const std::string &name,
                                     GLuint value_0, GLuint value_1,
                                     GLuint value_2, GLuint value_3) const;

        GLboolean setUniformValue1f(const std::string &name, GLfloat value) const;
        GLboolean setUniformValue2f(const std::string &name,
                                    GLfloat value_0, GLfloat value_1) const;
        GLboolean setUniformValue3f(const std::string &name,
                                    GLfloat value_0, GLfloat value_1,
                                    GLfloat value_2) const;
        GLboolean setUniformValue4f(const std::string &name, //(*@\label{src:ShaderProgram:setUniformValue4f:declaration}@*)
                                    GLfloat value_0, GLfloat value_1,
                                    GLfloat value_2, GLfloat value_3) const;

        GLboolean setUniformValue1iv(const std::string &name,
                                     GLsizei count, const GLint *value) const;
        GLboolean setUniformValue2iv(const std::string &name,
                                     GLsizei count, const GLint *value) const;
        GLboolean setUniformValue3iv(const std::string &name,
                                     GLsizei count, const GLint *value) const;
        GLboolean setUniformValue4iv(const std::string &name,
                                     GLsizei count, const GLint *value) const;

        GLboolean setUniformValue1uiv(const std::string &name,
                                      GLsizei count, const GLuint *value) const;
        GLboolean setUniformValue2uiv(const std::string &name,
                                      GLsizei count, const GLuint *value) const;
        GLboolean setUniformValue3uiv(const std::string &name,
                                      GLsizei count, const GLuint *value) const;
        GLboolean setUniformValue4uiv(const std::string &name,
                                      GLsizei count, const GLuint *value) const;

        GLboolean setUniformValue1fv(const std::string &name,
                                     GLsizei count, const GLfloat *value) const;
        GLboolean setUniformValue2fv(const std::string &name,
                                     GLsizei count, const GLfloat *value) const;
        GLboolean setUniformValue3fv(const std::string &name,
                                     GLsizei count, const GLfloat *value) const;
        GLboolean setUniformValue4fv(const std::string &name, //(*@\label{src:ShaderProgram:setUniformValue4fv:declaration}@*)
                                     GLsizei count, const GLfloat *value) const;

        GLboolean setUniformMatrix2fv(const std::string &name,
                                      GLsizei count, GLboolean transpose,
                                      const GLfloat *values) const;
        GLboolean setUniformMatrix3fv(const std::string &name,
                                      GLsizei count, GLboolean transpose,
                                      const GLfloat *values) const;
        GLboolean setUniformMatrix4fv(const std::string &name, //(*@\label{src:ShaderProgram:setUniformMatrix4fv:declaration}@*)
                                      GLsizei count, GLboolean transpose,
                                      const GLfloat *values) const;

        GLboolean setUniformMatrix2x3fv(const std::string &name,
                                        GLsizei count, GLboolean transpose,
                                        const GLfloat *values) const;
        GLboolean setUniformMatrix3x2fv(const std::string &name,
                                        GLsizei count, GLboolean transpose,
                                        const GLfloat *values) const;
        GLboolean setUniformMatrix2x4fv(const std::string &name,
                                        GLsizei count, GLboolean transpose,
                                        const GLfloat *values) const;
        GLboolean setUniformMatrix4x2fv(const std::string &name,
                                        GLsizei count, GLboolean transpose,
                                        const GLfloat *values) const;
        GLboolean setUniformMatrix3x4fv(const std::string &name,
                                        GLsizei count, GLboolean transpose,
                                        const GLfloat *values) const;
        GLboolean setUniformMatrix4x3fv(const std::string &name,
                                        GLsizei count, GLboolean transpose,
                                        const GLfloat *values) const;

        GLboolean setUniformColor(const std::string &name, const Color4 &color) const; //(*@\label{src:ShaderProgram:setUniformColor:declaration}@*)

        //(*@\Green{// The shader program assumes that uniform light source variables are of the type}@*)
        //(*@\Green{//}@*)
        //    struct LightSource
        //    {
        //        bool    enabled;
        //        vec4    position;
        //        vec4    half_vector;
        //        vec4    ambient;
        //        vec4    diffuse;
        //        vec4    specular;
        //        float   spot_cos_cutoff;
        //        float   constant_attenuation;
        //        float   linear_attenuation;
        //        float   quadratic_attenuation;
        //        vec3    spot_direction;
        //        float   spot_exponent;
        //    };
        GLboolean setUniformDirectionalLight(const std::string &name,
                                             const DirectionalLight &light) const;
        GLboolean setUniformPointLight(const std::string &name,
                                       const PointLight &light) const;
        GLboolean setUniformSpotlight(const std::string &name,
                                      const Spotlight &light) const;

        //(*@\Green{// The shader program assumes that uniform material variables are of the type}@*)
        //(*@\Green{//}@*)
        //     struct Material
        //     {
        //          vec4    ambient;
        //          vec4    diffuse;
        //          vec4    specular;
        //          vec4    emission;
        //          float   shininess;
        //     };
        GLboolean setUniformMaterial(const std::string &name,
                                     const Material &material) const;

        //(*@\Green{// The ambient and diffuse reflection coefficients of the}@*)
        //(*@\Green{//}@*)
        //     struct Material
        //     {
        //          vec4    ambient;
        //          vec4    diffuse;
        //          vec4    specular;
        //          vec4    emission;
        //          float   shininess;
        //     };
        //(*@\Green{//}@*)
        //(*@\Green{// will track the specified ambient\_and\_diffuse\_reflection\_coefficients, the emission color components}@*)
        //(*@\Green{// will be set to the fixed values (0.0f, 0.0f, 0.0f, 0.0f), while the specular reflection }@*)
        //(*@\Green{// coefficients and the shininess will correspond to the remaining user-specified values.}@*)
        GLboolean setUniformColorMaterial(
                const std::string &name,
                const Color4 &ambient_and_diffuse_reflection_coefficients,
                const Color4 &specular_reflection_coefficients =
                             Color4(0.3f, 0.3f, 0.3f, 0.8f),
                GLfloat shininess = 128.0f) const; //(*@\label{src:ShaderProgram:setUniformVariables:declaration:end}@*)

        //(*@\Green{// destructor}@*)
        virtual   ~ShaderProgram();
    };
}

#endif //(*@\Green{// SHADERPROGRAMS\_H}@*)
