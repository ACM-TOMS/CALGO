//----------------------------------------------------------------------------------
// File:        Core/Geometry/Surfaces/TriangleMeshes3.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef TRIANGLEMESHES3_H
#define TRIANGLEMESHES3_H

#include <GL/glew.h>

#include "../Coordinates/Cartesians3.h"
#include "../Coordinates/Colors4.h"
#include "../Coordinates/TCoordinates4.h"
#include "../../Shaders/ShaderPrograms.h"
#include "TriangularFaces.h"

#include <iostream>
#include <string>
#include <vector>

namespace cagd
{
    class TriangleMesh3
    {
        //(*@\Green{// a friend class that will be discussed later in Listings \mref{src:TensorProductSurfaces3.h} and \mref{src:TensorProductSurfaces3.cpp}}@*)
        friend class TensorProductSurface3;

        //(*@\Green{// overloaded friend input/output from/to stream operators}@*)
        friend std::ostream& operator <<(std::ostream& lhs, const TriangleMesh3& rhs);
        friend std::istream& operator >>(std::istream& lhs, TriangleMesh3& rhs);

    protected:
        //(*@\Green{// vertex buffer object identifiers and their common usage flag}@*)
        GLenum                      _usage_flag;
        GLuint                      _vbo_positions;
        GLuint                      _vbo_normals;
        GLuint                      _vbo_colors;
        GLuint                      _vbo_tex_coordinates;
        GLuint                      _vbo_indices;

        //(*@\Green{// corners of the bounding box}@*)
        Cartesian3                  _leftmost_position;
        Cartesian3                  _rightmost_position;

        //(*@\Green{// geometry}@*)
        std::vector<Cartesian3>     _position;
        std::vector<Cartesian3>     _normal;
        std::vector<Color4>         _color;
        std::vector<TCoordinate4>   _tex;
        std::vector<TriangularFace> _face;

    public:
        //(*@\Green{// default/special constructor}@*)
        TriangleMesh3(GLint vertex_count = 0, GLint face_count = 0,
                      GLenum usage_flag = GL_STATIC_DRAW);

        //(*@\Green{// copy constructor}@*)
        TriangleMesh3(const TriangleMesh3& mesh);

        //(*@\Green{// assignment operator}@*)
        TriangleMesh3& operator =(const TriangleMesh3& rhs);

        //(*@\Green{// Deletes all vertex buffer objects.}@*)
        GLvoid deleteVertexBufferObjects();

        //(*@\Green{// Geometry rendering methods.}@*)

        //(*@\Green{// The next rendering method will return GL\_FALSE if:\label{src:TrinagleMesh3:render:start}}@*)
        //(*@\Green{// \ \ - the given shader program is not active;}@*)
        //(*@\Green{// \ \ - the user-defined position, normal and texture attribute names cannot all be found in the}@*)
        //(*@\Green{// \ \ \ \ list of active attributes of the provided shader program, or they exist but all of}@*)
        //(*@\Green{// \ \ \ \ them have incorrect types;}@*)
        //(*@\Green{// \ \ - the render\_mode is incorrect (one should use either GL\_TRIANGLES or GL\_POINTS).}@*)
        GLboolean render(const ShaderProgram &program,
                         GLenum render_mode = GL_TRIANGLES) const;

        //(*@\Green{// If during rendering one intends to use shader program objects that are not instances of}@*)
        //(*@\Green{// our class ShaderProgram, then one should specify attribute locations associated with}@*)
        //(*@\Green{// positions of type vec3, unit normals of type vec3 and (possibly projective) texture}@*)
        //(*@\Green{// coordinates of type vec4.}@*)
        //(*@\Green{//}@*)
        //(*@\Green{// The next rendering method will return GL\_FALSE if:}@*)
        //(*@\Green{// \ \ - there is no active shader program object;}@*)
        //(*@\Green{// \ \ - the given position, normal and texture locations either cannot all be found in the list}@*)
        //(*@\Green{// \ \ \ \ of active attribute locations, or they exists but all of them have incorrect types;}@*)
        //(*@\Green{// \ \ - the render\_mode is incorrect (one should use either GL\_TRIANGLES or GL\_POINTS).}@*)
        GLboolean render(GLenum render_mode            = GL_TRIANGLES,
                         GLint  vec3_position_location = 0,
                         GLint  vec3_normal_location   = 1,
                         GLint  vec4_color_location    = 2,
                         GLint  vec4_texture_location  = 3) const; //(*@\label{src:TrinagleMesh3:render:end}@*)

        //(*@\Green{// Updates all vertex buffer objects.}@*)
        GLboolean updateVertexBufferObjects(GLenum usage_flag = GL_STATIC_DRAW);

        //(*@\Green{// Loads the geometry (i.e., the array of vertices and faces) stored in an OFF file,}@*)
        //(*@\Green{// by calculating at the same time the unit normal vectors associated with vertices.}@*)
        GLboolean loadFromOFF(const std::string &file_name,
                              GLboolean translate_and_scale_to_unit_cube = GL_FALSE);

        //(*@\Green{// for mapping vertex buffer objects}@*)
        GLfloat* mapPositionBuffer(GLenum access_flag = GL_READ_ONLY) const;
        GLfloat* mapNormalBuffer(GLenum access_flag = GL_READ_ONLY) const;
        GLfloat* mapColorBuffer(GLenum access_flag = GL_READ_ONLY) const;
        GLfloat* mapTextureBuffer(GLenum access_flag = GL_READ_ONLY) const;

        //(*@\Green{// for unmapping vertex buffer objects}@*)
        GLvoid unmapPositionBuffer() const;
        GLvoid unmapNormalBuffer() const;
        GLvoid unmapColorBuffer() const;
        GLvoid unmapTextureBuffer() const;

        //(*@\Green{// get properties of the geometry}@*)
        GLint vertexCount() const;
        GLint faceCount() const;

        Cartesian3& position(GLint index);
        const Cartesian3& position(GLint index) const;

        Cartesian3& normal(GLint index);
        const Cartesian3& normal(GLint index) const;

        Color4& color(GLint index);
        const Color4& color(GLint index) const;

        TCoordinate4& texture(GLint index);
        const TCoordinate4& texture(GLint index) const;

        TriangularFace& face(GLint index);
        const TriangularFace& face(GLint index) const;

        //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        virtual TriangleMesh3* clone() const;

        //(*@\Green{// destructor}@*)
        virtual ~TriangleMesh3();
    };
}

#endif //(*@\Green{// TRIANGLEMESHES3\_H}@*)
