//----------------------------------------------------------------------------------
// File:        Core/Geometry/Surfaces/TriangleMeshes3.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "TriangleMeshes3.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <limits>

using namespace std;

namespace cagd
{
    (*@\Green{// default/special constructor}@*)
    TriangleMesh3::TriangleMesh3(GLint vertex_count, GLint face_count, GLenum usage_flag):
        _usage_flag(usage_flag),
        _vbo_positions(0),
        _vbo_normals(0),
        _vbo_colors(0),
        _vbo_tex_coordinates(0),
        _vbo_indices(0),
        _position(vertex_count < 0 ? 0 : vertex_count),
        _normal(vertex_count < 0 ? 0 : vertex_count),
        _color(vertex_count < 0 ? 0 : vertex_count),
        _tex(vertex_count < 0 ? 0 : vertex_count),
        _face(face_count < 0 ? 0 : face_count)
    {
        assert("The number of vertices should be non-negative!" && vertex_count >= 0);

        assert("The number of triangular faces should be non-negative!" &&
               face_count >= 0);

        assert("Invalid usage flag!" &&
               (usage_flag == GL_STREAM_DRAW  || usage_flag == GL_STREAM_READ  ||
                usage_flag == GL_STREAM_COPY  ||
                usage_flag == GL_DYNAMIC_DRAW || usage_flag == GL_DYNAMIC_READ ||
                usage_flag == GL_DYNAMIC_COPY ||
                usage_flag == GL_STATIC_DRAW  || usage_flag == GL_STATIC_READ  ||
                usage_flag == GL_STATIC_COPY));
    }

    (*@\Green{// copy constructor}@*)
    TriangleMesh3::TriangleMesh3(const TriangleMesh3 &mesh):
            _usage_flag(mesh._usage_flag),
            _vbo_positions(0),
            _vbo_normals(0),
            _vbo_colors(0),
            _vbo_tex_coordinates(0),
            _vbo_indices(0),
            _position(mesh._position),
            _normal(mesh._normal),
            _color(mesh._color),
            _tex(mesh._tex),
            _face(mesh._face)
    {
        if (mesh._vbo_positions &&
            mesh._vbo_normals &&
            mesh._vbo_colors &&
            mesh._vbo_tex_coordinates &&
            mesh._vbo_indices)
        {
            updateVertexBufferObjects(mesh._usage_flag);
        }
    }

    (*@\Green{// assignment operator}@*)
    TriangleMesh3& TriangleMesh3::operator =(const TriangleMesh3 &rhs)
    {
        if (this != &rhs)
        {
            deleteVertexBufferObjects();

            _usage_flag = rhs._usage_flag;
            _position   = rhs._position;
            _normal     = rhs._normal;
            _color      = rhs._color;
            _tex        = rhs._tex;
            _face       = rhs._face;

            if (rhs._vbo_positions &&
                rhs._vbo_normals &&
                rhs._vbo_colors &&
                rhs._vbo_tex_coordinates &&
                rhs._vbo_indices)
            {
                updateVertexBufferObjects(_usage_flag);
            }
        }

        return *this;
    }

    (*@\Green{// Deletes all vertex buffer objects.}@*)
    GLvoid TriangleMesh3::deleteVertexBufferObjects()
    {
        if (_vbo_positions)
        {
            glDeleteBuffers(1, &_vbo_positions);
            _vbo_positions = 0;
        }

        if (_vbo_normals)
        {
            glDeleteBuffers(1, &_vbo_normals);
            _vbo_normals = 0;
        }

        if (_vbo_colors)
        {
            glDeleteBuffers(1, &_vbo_colors);
            _vbo_colors = 0;
        }

        if (_vbo_tex_coordinates)
        {
            glDeleteBuffers(1, &_vbo_tex_coordinates);
            _vbo_tex_coordinates = 0;
        }

        if (_vbo_indices)
        {
            glDeleteBuffers(1, &_vbo_indices);
            _vbo_indices = 0;
        }
    }

    (*@\Green{// Geometry rendering methods.}@*)

    (*@\Green{// The next rendering method will return GL\_FALSE if:}@*)
    (*@\Green{// \ \ - the given shader program is not active;}@*)
    (*@\Green{// \ \ - the user-defined position, normal and texture attribute names cannot all be found in the}@*)
    (*@\Green{// \ \ \ \ list of active attributes of the provided shader program, or they exist but all of}@*)
    (*@\Green{// \ \ \ \ them have incorrect types;}@*)
    (*@\Green{// \ \ - the render\_mode is incorrect (one should use either GL\_TRIANGLES or GL\_POINTS).}@*)
    GLboolean TriangleMesh3::render(const ShaderProgram &program,
                                    GLenum render_mode) const
    {
        GLint current_program;
        glGetIntegerv(GL_CURRENT_PROGRAM, &current_program);

        if (current_program != (GLint)program._id)
        {
            return GL_FALSE;
        }

        GLint position = program.positionAttributeLocation();
        GLint normal   = program.normalAttributeLocation();
        GLint color    = program.colorAttributeLocation();
        GLint texture  = program.textureAttributeLocation();

        if (!_vbo_positions || !_vbo_normals || !_vbo_colors || !_vbo_tex_coordinates ||
            !_vbo_indices || (position < 0 && normal < 0 && color < 0 && texture < 0))
        {
            return GL_FALSE;
        }

        if (render_mode != GL_TRIANGLES && render_mode != GL_POINTS)
        {
            return GL_FALSE;
        }

        (*@\Green{// enable array of texture, color, normal and position attributes}@*)
        if (texture >= 0)
        {
            glEnableVertexAttribArray(texture);
            (*@\Green{// activate the VBO of texture coordinates}@*)
            glBindBuffer(GL_ARRAY_BUFFER, _vbo_tex_coordinates);
            (*@\Green{// specify the location and data format of texture coordinates}@*)
            glVertexAttribPointer(texture, 4, GL_FLOAT, GL_FALSE, 0, nullptr);
        }

        if (color >= 0)
        {
            glEnableVertexAttribArray(color);
            (*@\Green{// activate the VBO of colors}@*)
            glBindBuffer(GL_ARRAY_BUFFER, _vbo_colors);
            (*@\Green{// specify the location and data format of colors}@*)
            glVertexAttribPointer(color, 4, GL_FLOAT, GL_FALSE, 0, nullptr);
        }

        if (normal >= 0)
        {
            glEnableVertexAttribArray(normal);
            (*@\Green{// activate the VBO of normal vectors}@*)
            glBindBuffer(GL_ARRAY_BUFFER, _vbo_normals);
            (*@\Green{// specify the location and data format of normal vectors}@*)
            glVertexAttribPointer(normal, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        }

        if (position >= 0)
        {
            glEnableVertexAttribArray(position);
            (*@\Green{// activate the VBO of vertices}@*)
            glBindBuffer(GL_ARRAY_BUFFER, _vbo_positions);
            (*@\Green{// specify the location and data format of vertices}@*)
            glVertexAttribPointer(position, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        }

        (*@\Green{// activate the element array buffer for indexed vertices of triangular faces}@*)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _vbo_indices);

        (*@\Green{// render primitives}@*)
        glDrawElements(render_mode, 3 * (GLsizei)_face.size(), GL_UNSIGNED_INT, nullptr);

        (*@\Green{// disable individual attribute arrays}@*)
        if (position >= 0)
        {
            glDisableVertexAttribArray(position);
            glVertexAttrib3f(position, 0.0f, 0.0f, 0.0f);
        }

        if (normal >= 0)
        {
            glDisableVertexAttribArray(normal);
            glVertexAttrib3f(normal, 0.0f, 0.0f, 0.0f);
        }

        if (color >= 0)
        {
            glDisableVertexAttribArray(color);
            glVertexAttrib4f(color, 0.0f, 0.0f, 0.0f, 1.0f);
        }

        if (texture >= 0)
        {
            glDisableVertexAttribArray(texture);
            glVertexAttrib4f(texture, 0.0f, 0.0f, 0.0f, 1.0f);
        }

        (*@\Green{// unbind any buffer object previously bound and restore client memory usage}@*)
        (*@\Green{// for these buffer object targets}@*)
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

        return GL_TRUE;
    }

    (*@\Green{// If during rendering one intends to use shader program objects that are not instances of}@*)
    (*@\Green{// our class ShaderProgram, then one should specify attribute locations associated with}@*)
    (*@\Green{// positions of type vec3, unit normals of type vec3 and (possibly projective) texture}@*)
    (*@\Green{// coordinates of type vec4.}@*)
    (*@\Green{//}@*)
    (*@\Green{// The next rendering method will return GL\_FALSE if:}@*)
    (*@\Green{// \ \ - there is no active shader program object;}@*)
    (*@\Green{// \ \ - the given position, normal and texture locations either cannot all be found in the list}@*)
    (*@\Green{// \ \ \ \ of active attribute locations, or they exists but all of them have incorrect types;}@*)
    (*@\Green{// \ \ - the render\_mode is incorrect (one should use either GL\_TRIANGLES or GL\_POINTS).}@*)
    GLboolean TriangleMesh3::render(
            GLenum render_mode,
            GLint vec3_position_location,
            GLint vec3_normal_location,
            GLint vec4_color_location,
            GLint vec4_texture_location) const
    {
        if (!_vbo_positions || !_vbo_normals || !_vbo_colors || !_vbo_tex_coordinates ||
            !_vbo_indices ||
            (vec3_position_location < 0 && vec3_normal_location < 0 &&
             vec4_color_location < 0 && vec4_texture_location < 0))
        {
            return GL_FALSE;
        }

        if (render_mode != GL_TRIANGLES && render_mode != GL_POINTS)
        {
            return GL_FALSE;
        }

        GLint current_program = 0;
        glGetIntegerv(GL_CURRENT_PROGRAM, &current_program);

        if (!current_program)
        {
            return GL_FALSE;
        }

        GLint attribute_count;
        glGetProgramiv(current_program, GL_ACTIVE_ATTRIBUTES, &attribute_count);

        if (!attribute_count)
        {
            return GL_FALSE;
        }

        GLsizei max_attribute_name_length;
        glGetProgramiv(current_program, GL_ACTIVE_ATTRIBUTE_MAX_LENGTH,
                       &max_attribute_name_length);

        GLchar *attribute_name_data = new GLchar[max_attribute_name_length];

        GLint     given_attribute_count = 4;

        GLint     given_locations[] = {vec3_position_location, vec3_normal_location,
                                       vec4_color_location, vec4_texture_location};

        GLenum    expected_types[]  = {GL_FLOAT_VEC3, GL_FLOAT_VEC3,
                                       GL_FLOAT_VEC4, GL_FLOAT_VEC4};

        GLboolean given_locations_exist_and_have_proper_types[] = {GL_FALSE, GL_FALSE,
                                                                   GL_FALSE, GL_FALSE};

        for (GLint attribute = 0;
             attribute < attribute_count &&
             (!given_locations_exist_and_have_proper_types[0] ||
              !given_locations_exist_and_have_proper_types[1] ||
              !given_locations_exist_and_have_proper_types[2] ||
              !given_locations_exist_and_have_proper_types[3]);
             attribute++)
        {
            GLsizei actual_length = 0;
            GLint   array_size    = 0;
            GLenum  type          = 0;

            glGetActiveAttrib(current_program, attribute, max_attribute_name_length,
                              &actual_length, &array_size, &type, attribute_name_data);
            string name(&attribute_name_data[0], actual_length);
            GLint  location = glGetAttribLocation(current_program, name.c_str());

            for (GLint i = 0; i < given_attribute_count; i++)
            {
                if (type == expected_types[i] && given_locations[i] == location)
                {
                    given_locations_exist_and_have_proper_types[i] = GL_TRUE;
                }
            }
        }

        delete[] attribute_name_data;

        if (!given_locations_exist_and_have_proper_types[0] &&
            !given_locations_exist_and_have_proper_types[1] &&
            !given_locations_exist_and_have_proper_types[2] &&
            !given_locations_exist_and_have_proper_types[3])
        {
            return GL_FALSE;
        }

        (*@\Green{// enable array of position, normal, color and texture attributes}@*)
        for (GLint i = 0; i < given_attribute_count; i++)
        {
            if (given_locations_exist_and_have_proper_types[i] &&
                (given_locations[i] >= 0))
            {
                glEnableVertexAttribArray(given_locations[i]);
            }
        }

        GLuint vbos[]             = {_vbo_positions, _vbo_normals,
                                     _vbo_colors, _vbo_tex_coordinates};
        GLint  component_counts[] = {3, 3, 4, 4};

        for (GLint i = 0; i < given_attribute_count; i++)
        {
            if (given_locations_exist_and_have_proper_types[i] &&
                (given_locations[i] >= 0))
            {
                (*@\Green{// activate the $i$th VBO}@*)
                glBindBuffer(GL_ARRAY_BUFFER, vbos[i]);

                (*@\Green{// and specify its location and data format}@*)
                glVertexAttribPointer(given_locations[i], component_counts[i],
                                      GL_FLOAT, GL_FALSE, 0, nullptr);
            }
        }

        (*@\Green{// activate the element array buffer for indexed vertices of triangular faces}@*)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _vbo_indices);

        (*@\Green{// render primitives}@*)
        glDrawElements(render_mode, 3 * (GLsizei)_face.size(), GL_UNSIGNED_INT, nullptr);

        (*@\Green{// unbind any buffer object previously bound and restore client memory usage}@*)
        (*@\Green{// for these buffer object targets}@*)
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

        (*@\Green{// disable individual attribute arrays}@*)
        for (GLint i = 0; i < given_attribute_count; i++)
        {
            if (given_locations_exist_and_have_proper_types[i] &&
                (given_locations[i] >= 0))
            {
                glDisableVertexAttribArray(given_locations[i]);

                switch (expected_types[i])
                {
                case GL_FLOAT_VEC3:
                    glVertexAttrib3f(given_locations[i], 0.0f, 0.0f, 0.0f);
                    break;

                default:
                    glVertexAttrib4f(given_locations[i], 0.0f, 0.0f, 0.0f, 1.0f);
                    break;
                }
            }
        }

        return GL_TRUE;
    }

    (*@\Green{// Updates all vertex buffer objects.}@*)
    GLboolean TriangleMesh3::updateVertexBufferObjects(GLenum usage_flag)
    {
        if (usage_flag != GL_STREAM_DRAW  && usage_flag != GL_STREAM_READ  &&
            usage_flag != GL_STREAM_COPY  &&
            usage_flag != GL_STATIC_DRAW  && usage_flag != GL_STATIC_READ  &&
            usage_flag != GL_STATIC_COPY  &&
            usage_flag != GL_DYNAMIC_DRAW && usage_flag != GL_DYNAMIC_READ &&
            usage_flag != GL_DYNAMIC_COPY)
        {
            return GL_FALSE;
        }

        (*@\Green{// updating the usage flag}@*)
        _usage_flag = usage_flag;

        (*@\Green{// deleting old vertex buffer objects}@*)
        deleteVertexBufferObjects();

        (*@\Green{// creating vertex buffer objects of mesh vertices, unit normal vectors, colors, texture coordinates}@*)
        (*@\Green{// and element indices}@*)
        glGenBuffers(1, &_vbo_positions);

        if (!_vbo_positions)
        {
            return GL_FALSE;
        }

        glGenBuffers(1, &_vbo_normals);

        if (!_vbo_normals)
        {
            glDeleteBuffers(1, &_vbo_positions);
            _vbo_positions = 0;

            return GL_FALSE;
        }

        glGenBuffers(1, &_vbo_colors);

        if (!_vbo_colors)
        {
            glDeleteBuffers(1, &_vbo_positions);
            _vbo_positions = 0;

            glDeleteBuffers(1, &_vbo_normals);
            _vbo_normals = 0;

            return GL_FALSE;
        }

        glGenBuffers(1, &_vbo_tex_coordinates);
        if (!_vbo_tex_coordinates)
        {
            glDeleteBuffers(1, &_vbo_positions);
            _vbo_positions = 0;

            glDeleteBuffers(1, &_vbo_normals);
            _vbo_normals = 0;

            glDeleteBuffers(1, &_vbo_colors);
            _vbo_colors = 0;

            return GL_FALSE;
        }

        glGenBuffers(1, &_vbo_indices);
        if (!_vbo_indices)
        {
            glDeleteBuffers(1, &_vbo_positions);
            _vbo_positions = 0;

            glDeleteBuffers(1, &_vbo_normals);
            _vbo_normals = 0;

            glDeleteBuffers(1, &_vbo_colors);
            _vbo_colors = 0;

            glDeleteBuffers(1, &_vbo_tex_coordinates);
            _vbo_tex_coordinates = 0;

            return GL_FALSE;
        }

        (*@\Green{// For efficiency reasons we convert all GLdouble coordinates to GLfloat ones: we will use }@*)
        (*@\Green{// auxiliar pointers for buffer data loading and functions glMapBuffer/glUnmapBuffer.}@*)
        (*@\Green{// Note that, multiple buffers can be mapped simultaneously.}@*)

        GLsizeiptr vertex_byte_size = 3 * _position.size() * sizeof(GLfloat);

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_positions);
        glBufferData(GL_ARRAY_BUFFER, vertex_byte_size, nullptr, _usage_flag);

        GLfloat *vertex_coordinate = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER,
                                                           GL_WRITE_ONLY);

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_normals);
        glBufferData(GL_ARRAY_BUFFER, vertex_byte_size, nullptr, _usage_flag);

        GLfloat *normal_coordinate = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER,
                                                           GL_WRITE_ONLY);

        for (vector<Cartesian3>::const_iterator
             vit = _position.begin(),
             nit = _normal.begin(); vit != _position.end(); vit++, nit++)
        {
            for (GLint component = 0; component < 3; component++,
                 vertex_coordinate++, normal_coordinate++)
            {
                *vertex_coordinate = (GLfloat)(*vit)[component];
                *normal_coordinate = (GLfloat)(*nit)[component];
            }
        }

        GLsizeiptr color_or_tex_byte_size = 4 * _color.size() * sizeof(GLfloat);

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_colors);
        glBufferData(GL_ARRAY_BUFFER, color_or_tex_byte_size, nullptr, _usage_flag);
        GLfloat *color_coordinate = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

        memcpy(color_coordinate, &_color[0][0], color_or_tex_byte_size);

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_tex_coordinates);
        glBufferData(GL_ARRAY_BUFFER, color_or_tex_byte_size, nullptr, _usage_flag);
        GLfloat *tex_coordinate = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

        memcpy(tex_coordinate, &_tex[0][0], color_or_tex_byte_size);

        GLsizeiptr index_byte_size = 3 * (GLsizeiptr)(_face.size() * sizeof(GLint));

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _vbo_indices);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, index_byte_size, nullptr, _usage_flag);
        GLuint *element = (GLuint*)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);

        for (vector<TriangularFace>::const_iterator fit = _face.begin();
             fit != _face.end(); fit++)
        {
            for (GLint node = 0; node < 3; node++, element++)
            {
                *element = (*fit)[node];
            }
        }

        (*@\Green{// unmap all VBOs}@*)
        glBindBuffer(GL_ARRAY_BUFFER, _vbo_positions);
        if (!glUnmapBuffer(GL_ARRAY_BUFFER))
            return GL_FALSE;

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_normals);
        if (!glUnmapBuffer(GL_ARRAY_BUFFER))
            return GL_FALSE;

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_colors);
        if (!glUnmapBuffer(GL_ARRAY_BUFFER))
            return GL_FALSE;

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_tex_coordinates);
        if (!glUnmapBuffer(GL_ARRAY_BUFFER))
            return GL_FALSE;

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _vbo_indices);
        if (!glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER))
            return GL_FALSE;

        (*@\Green{// unbind any buffer object previously bound and restore client memory usage for these buffer object}@*)
        (*@\Green{// targets}@*)
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

        return GL_TRUE;
    }

    (*@\Green{// loads the geometry (i.e., the array of vertices and faces) stored in an OFF file}@*)
    (*@\Green{// at the same time calculates the unit normal vectors associated with vertices}@*)
    GLboolean TriangleMesh3::loadFromOFF(
            const string &file_name, GLboolean translate_and_scale_to_unit_cube)
    {
        fstream f(file_name.c_str(), ios_base::in);

        if (!f || !f.good())
        {
            return GL_FALSE;
        }

        (*@\Green{// loading the header}@*)
        string header;

        f >> header;

        if (header != "OFF")
        {
            return GL_FALSE;
        }

        (*@\Green{// loading number of vertices, faces, and edges}@*)
        GLint vertex_count, face_count, edge_count;

        f >> vertex_count >> face_count >> edge_count;

        (*@\Green{// allocating memory for vertices, unit normal vectors, colors, texture coordinates and faces}@*)
        _position.resize(vertex_count);
        _normal.resize(vertex_count);
        _color.resize(vertex_count);
        _tex.resize(vertex_count);
        _face.resize(face_count);

        (*@\Green{// initializing the leftmost and rightmost corners of the bounding box}@*)
        _leftmost_position.x() = _leftmost_position.y() = _leftmost_position.z() =
                 numeric_limits<GLdouble>::max();
        _rightmost_position.x() = _rightmost_position.y() = _rightmost_position.z() =
                -numeric_limits<GLdouble>::max();

        (*@\Green{// loading vertices and correcting the leftmost and rightmost corners of the bounding box}@*)
        for (vector<Cartesian3>::iterator vit = _position.begin();
             vit != _position.end(); vit++)
        {
            f >> *vit;

            if (vit->x() < _leftmost_position.x()) {_leftmost_position.x() = vit->x();}
            if (vit->y() < _leftmost_position.y()) {_leftmost_position.y() = vit->y();}
            if (vit->z() < _leftmost_position.z()) {_leftmost_position.z() = vit->z();}

            if (vit->x() > _rightmost_position.x()) {_rightmost_position.x() = vit->x();}
            if (vit->y() > _rightmost_position.y()) {_rightmost_position.y() = vit->y();}
            if (vit->z() > _rightmost_position.z()) {_rightmost_position.z() = vit->z();}
        }

        (*@\Green{// if we do not want to preserve the original positions and coordinates of vertices}@*)
        if (translate_and_scale_to_unit_cube)
        {
            Cartesian3 diagonal(_rightmost_position);
            diagonal -= _leftmost_position;

            GLdouble scale = 1.0 / max(diagonal[0], max(diagonal[1], diagonal[2]));

            Cartesian3 middle(_leftmost_position);
            middle += _rightmost_position;
            middle *= 0.5;
            for (vector<Cartesian3>::iterator vit = _position.begin();
                 vit != _position.end(); vit++)
            {
                *vit -= middle;
                *vit *= scale;
            }
        }

        (*@\Green{// loading faces}@*)
        for (vector<TriangularFace>::iterator fit = _face.begin();
             fit != _face.end(); fit++)
        {
            f >> *fit;
        }

        (*@\Green{// calculating average unit normal vectors associated with vertices}@*)
        for (vector<TriangularFace>::const_iterator fit = _face.begin();
             fit != _face.end(); fit++)
        {
            Cartesian3 n = _position[(*fit)[1]];
            n -= _position[(*fit)[0]];

            Cartesian3 p = _position[(*fit)[2]];
            p -= _position[(*fit)[0]];

            n ^= p;

            for (GLint node = 0; node < 3; node++)
            {
                _normal[(*fit)[node]] += n;
            }
        }

        for (vector<Cartesian3>::iterator nit = _normal.begin();
             nit != _normal.end(); nit++)
        {
            nit->normalize();
        }

        f.close();

        return GL_TRUE;
    }

    (*@\Green{// for mapping vertex buffer objects}@*)
    GLfloat* TriangleMesh3::mapPositionBuffer(GLenum access_flag) const
    {
        if (access_flag != GL_READ_ONLY && access_flag != GL_WRITE_ONLY &&
            access_flag != GL_READ_WRITE)
        {
            return nullptr;
        }

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_positions);
        GLfloat* result = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, access_flag);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        return result;
    }

    GLfloat* TriangleMesh3::mapNormalBuffer(GLenum access_flag) const
    {
        if (access_flag != GL_READ_ONLY && access_flag != GL_WRITE_ONLY &&
            access_flag != GL_READ_WRITE)
        {
            return nullptr;
        }

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_normals);
        GLfloat* result = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, access_flag);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        return result;
    }

    GLfloat* TriangleMesh3::mapColorBuffer(GLenum access_flag) const
    {
        if (access_flag != GL_READ_ONLY && access_flag != GL_WRITE_ONLY &&
            access_flag != GL_READ_WRITE)
        {
            return nullptr;
        }

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_colors);
        GLfloat* result = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, access_flag);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        return result;
    }

    GLfloat* TriangleMesh3::mapTextureBuffer(GLenum access_flag) const
    {
        if (access_flag != GL_READ_ONLY && access_flag != GL_WRITE_ONLY &&
            access_flag != GL_READ_WRITE)
        {
            return nullptr;
        }

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_tex_coordinates);
        GLfloat* result = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, access_flag);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        return result;
    }

    (*@\Green{// for unmapping vertex buffer objects}@*)
    GLvoid TriangleMesh3::unmapPositionBuffer() const
    {
        glBindBuffer(GL_ARRAY_BUFFER, _vbo_positions);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    GLvoid TriangleMesh3::unmapNormalBuffer() const
    {
        glBindBuffer(GL_ARRAY_BUFFER, _vbo_normals);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    GLvoid TriangleMesh3::unmapColorBuffer() const
    {
        glBindBuffer(GL_ARRAY_BUFFER, _vbo_normals);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    GLvoid TriangleMesh3::unmapTextureBuffer() const
    {
        glBindBuffer(GL_ARRAY_BUFFER, _vbo_tex_coordinates);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    (*@\Green{// get properties of the geometry}@*)
    GLint TriangleMesh3::vertexCount() const
    {
        return (GLint)_position.size();
    }

    GLint TriangleMesh3::faceCount() const
    {
        return (GLint)_face.size();
    }

    Cartesian3& TriangleMesh3::position(GLint index)
    {
        assert("The given position index is out of bounds!" &&
               (index >= 0 && index < (GLint)_position.size()));
        return _position[index];
    }

    const Cartesian3& TriangleMesh3::position(GLint index) const
    {
        assert("The given position index is out of bounds!" &&
               (index >= 0 && index < (GLint)_position.size()));
        return _position[index];
    }

    Cartesian3& TriangleMesh3::normal(GLint index)
    {
        assert("The given normal index is out of bounds!" &&
               (index >= 0 && index < (GLint)_normal.size()));
        return _normal[index];
    }

    const Cartesian3& TriangleMesh3::normal(GLint index) const
    {
        assert("The given normal index is out of bounds!" &&
               (index >= 0 && index < (GLint)_normal.size()));
        return _normal[index];
    }

    Color4& TriangleMesh3::color(GLint index)
    {
        assert("The given color index is out of bounds!" &&
               (index >= 0 && index < (GLint)_color.size()));
        return _color[index];
    }

    const Color4& TriangleMesh3::color(GLint index) const
    {
        assert("The given color index is out of bounds!" &&
               (index >= 0 && index < (GLint)_color.size()));
        return _color[index];
    }

    TCoordinate4& TriangleMesh3::texture(GLint index)
    {
        assert("The given texture index is out of bounds!" &&
               (index >= 0 && index < (GLint)_tex.size()));
        return _tex[index];
    }

    const TCoordinate4& TriangleMesh3::texture(GLint index) const
    {
        assert("The given texture index is out of bounds!" &&
               (index >= 0 && index < (GLint)_tex.size()));
        return _tex[index];
    }

    TriangularFace& TriangleMesh3::face(GLint index)
    {
        assert("The given face index is out of bounds!" &&
               (index >= 0 && index < (GLint)_face.size()));
        return _face[index];
    }

    const TriangularFace& TriangleMesh3::face(GLint index) const
    {
        assert("The given face index is out of bounds!" &&
               (index >= 0 && index < (GLint)_face.size()));
        return _face[index];
    }

    (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    TriangleMesh3* TriangleMesh3::clone() const
    {
        return new (nothrow) TriangleMesh3(*this);
    }

    (*@\Green{// destructor}@*)
    TriangleMesh3::~TriangleMesh3()
    {
        deleteVertexBufferObjects();
    }
}
