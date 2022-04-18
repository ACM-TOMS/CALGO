//----------------------------------------------------------------------------------
// File:        Core/Shaders/ShaderPrograms.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "ShaderPrograms.h"
#include "../Utilities.h"
#include <fstream>

using namespace std;

namespace cagd
{
    ShaderProgram::Shader::Shader(): _id(0), _compile_status(GL_FALSE), _type(VERTEX)
    {
    }

    ShaderProgram::Shader::Shader(const ShaderProgram::Shader &shader):
        _id(0), _compile_status(GL_FALSE), _type(shader._type), _source(shader._source)
    {
        switch (_type)
        {
        case VERTEX:
            _id = glCreateShader(GL_VERTEX_SHADER);
            break;

        case FRAGMENT:
            _id = glCreateShader(GL_FRAGMENT_SHADER);
            break;

        case COMPUTE:
            _id = glCreateShader(GL_COMPUTE_SHADER);
            break;

        case GEOMETRY:
            _id = glCreateShader(GL_GEOMETRY_SHADER);
            break;

        case TESS_CONTROL:
            _id = glCreateShader(GL_TESS_CONTROL_SHADER);
            break;

        case TESS_EVALUATION:
            _id = glCreateShader(GL_TESS_EVALUATION_SHADER);
            break;
        }

        const GLchar *pointer_to_source = _source.c_str();
        glShaderSource(_id, 1, &pointer_to_source, nullptr);

        if (shader._compile_status)
        {
            glCompileShader(_id);
            glGetShaderiv(_id, GL_COMPILE_STATUS, &_compile_status);
        }
        else
        {
            _compile_status = GL_FALSE;
        }
    }

    ShaderProgram::Shader& ShaderProgram::Shader::operator =(
            const ShaderProgram::Shader &rhs)
    {
        if (this != &rhs)
        {
            flagForDeletion();

            _type   = rhs._type;
            _source = rhs._source;

            switch (_type)
            {
            case VERTEX:
                _id = glCreateShader(GL_VERTEX_SHADER);
                break;

            case FRAGMENT:
                _id = glCreateShader(GL_FRAGMENT_SHADER);
                break;

            case COMPUTE:
                _id = glCreateShader(GL_COMPUTE_SHADER);
                break;

            case GEOMETRY:
                _id = glCreateShader(GL_GEOMETRY_SHADER);
                break;

            case TESS_CONTROL:
                _id = glCreateShader(GL_TESS_CONTROL_SHADER);
                break;

            case TESS_EVALUATION:
                _id = glCreateShader(GL_TESS_EVALUATION_SHADER);
                break;
            }

            const GLchar *pointer_to_source = _source.c_str();
            glShaderSource(_id, 1, &pointer_to_source, nullptr);

            if (rhs._compile_status)
            {
                glCompileShader(_id);
                glGetShaderiv(_id, GL_COMPILE_STATUS, &_compile_status);
            }
            else
            {
                _compile_status = GL_FALSE;
            }
        }

        return *this;
    }

    GLboolean ShaderProgram::Shader::createLoadAndCompileFromSourceFile(
            Type type, const std::string &file_name,
            GLboolean logging_is_enabled, ostream &output)
    {
        flagForDeletion();

        _type = type;

        switch (_type)
        {
        case VERTEX:
            _id = glCreateShader(GL_VERTEX_SHADER);
            break;

        case FRAGMENT:
            _id = glCreateShader(GL_FRAGMENT_SHADER);
            break;

        case COMPUTE:
            _id = glCreateShader(GL_COMPUTE_SHADER);
            break;

        case GEOMETRY:
            _id = glCreateShader(GL_GEOMETRY_SHADER);
            break;

        case TESS_CONTROL:
            _id = glCreateShader(GL_TESS_CONTROL_SHADER);
            break;

        case TESS_EVALUATION:
            _id = glCreateShader(GL_TESS_EVALUATION_SHADER);
            break;
        }

        string shader_type[6] = {"vertex", "fragment", "compute",
                                 "geometry", "tess-control", "tess-evaluation"};

        if (!_id)
        {
            if (logging_is_enabled)
            {
                output << "Could not create an empty " << shader_type[type]
                       << " shader object!" << endl;
            }
            return GL_FALSE;
        }

        fstream file(file_name.c_str(), ios_base::in);

        if (!file)
        {
            if (logging_is_enabled)
            {
                output << "The given file " << file_name << " does not exist!" << endl;
            }

            file.close();

            return GL_FALSE;
        }

        if (logging_is_enabled)
        {
            string caption = "Source of the " + shader_type[type] + " shader";
            string delimeter(caption.length(), '-');
            output << caption << endl << delimeter << endl;
        }

        _source.clear();

        string aux;
        while (!file.eof())
        {
            getline(file, aux, '\n');
            _source += aux + '\n';

            if (logging_is_enabled)
            {
                output << "\t" << aux << endl;
            }
        }

        file.close();

        if (logging_is_enabled)
        {
            output << endl;
        }

        const GLchar *pointer_to_source = _source.c_str();
        glShaderSource(_id, 1, &pointer_to_source, nullptr);

        glCompileShader(_id);
        glGetShaderiv(_id, GL_COMPILE_STATUS, &_compile_status);

        if (!_compile_status)
        {
            if (logging_is_enabled)
            {
                GLint info_log_length;
                glGetShaderiv(_id, GL_INFO_LOG_LENGTH, &info_log_length);

                if (info_log_length > 0)
                {
                    GLchar *info_log = new GLchar[info_log_length];
                    glGetShaderInfoLog(_id, info_log_length, nullptr, info_log);

                    output << "\t\\begin{" << shader_type[type]
                           << " shader information log}" << endl << "\t\tid = " << _id
                           << ", name = " << file_name << endl;
                    output <<  "\t\t" << info_log << endl;
                    output << "\t\\end{" << shader_type[type]
                           << " shader information log}" << endl << endl;

                    delete[] info_log;
                }
            }

            glDeleteShader(_id);
            return GL_FALSE;
        }

        if (logging_is_enabled)
        {
            output << "The given " << shader_type[type]
                   << " shader was successfully created, loaded and compiled."
                   << endl << endl;
        }

        return GL_TRUE;
    }

    GLvoid ShaderProgram::Shader::flagForDeletion() const
    {
        if (_id)
        {
            glDeleteShader(_id);
        }
    }

    ShaderProgram::Shader::~Shader()
    {
        flagForDeletion();
    }

    ShaderProgram::ActiveVariable::ActiveVariable(
            GLenum type, GLint array_size, GLint location):
        _type(type), _array_size(array_size), _location(location)
    {
    }

    GLenum ShaderProgram::ActiveVariable::type() const
    {
        return _type;
    }

    const GLchar *ShaderProgram::ActiveVariable::typeDescription() const
    {
        return openGLTypeToString(_type);
    }

    GLint ShaderProgram::ActiveVariable::arraySize() const
    {
        return _array_size;
    }

    GLint ShaderProgram::ActiveVariable::location() const
    {
        return _location;
    }

    ShaderProgram::ShaderProgram():
        _id(0), _link_status(GL_FALSE),
        _position_attribute_name("position"),
        _normal_attribute_name("normal"),
        _color_attribute_name("color"),
        _texture_attribute_name("texture")
    {        
    }

    ShaderProgram::ShaderProgram(const ShaderProgram &program):
        _id(0), _link_status(GL_FALSE),
        _shader(program._shader),
        _attribute(program._attribute),
        _uniform(program._uniform),
        _position_attribute_name(program._position_attribute_name),
        _normal_attribute_name(program._normal_attribute_name),
        _color_attribute_name(program._color_attribute_name),
        _texture_attribute_name(program._texture_attribute_name)
    {
        _id = glCreateProgram();

        if (_id)
        {
            for (GLuint i = 0; i < (GLuint)_shader.size(); i++)
            {
                if (_shader[i]._id)
                {
                    glAttachShader(_id, _shader[i]._id);
                }
            }

            if (program._link_status)
            {
                glLinkProgram(_id);
                glGetProgramiv(_id, GL_LINK_STATUS, &_link_status);
            }
        }
    }

    ShaderProgram& ShaderProgram::operator =(const ShaderProgram &rhs)
    {
        if (this != &rhs)
        {
            for (GLuint i = 0; i < _shader.size(); i++)
            {
                if (_shader[i]._id)
                {
                    glDetachShader(_id, _shader[i]._id);
                    _shader[i].flagForDeletion();
                }
            }

            if (_id)
            {
                glDeleteProgram(_id);
                _id = 0;
            }

            _link_status             = GL_FALSE;
            _shader                  = rhs._shader;
            _attribute               = rhs._attribute;
            _uniform                 = rhs._uniform;
            _position_attribute_name = rhs._position_attribute_name;
            _normal_attribute_name   = rhs._normal_attribute_name;
            _color_attribute_name    = rhs._color_attribute_name;
            _texture_attribute_name  = rhs._texture_attribute_name;

            _id = glCreateProgram();

            if (_id)
            {
                for (GLuint i = 0; i < (GLuint)_shader.size(); i++)
                {
                    if (_shader[i]._id)
                    {
                        glAttachShader(_id, _shader[i]._id);
                    }
                }

                if (rhs._link_status)
                {
                    glLinkProgram(_id);
                    glGetProgramiv(_id, GL_LINK_STATUS, &_link_status);
                }
            }
        }

        return *this;
    }

    GLboolean ShaderProgram::attachNewShaderFromSourceFile(
            Shader::Type type, const string &file_name,
            GLboolean logging_is_enabled, ostream &output)
    {
        if (!_id)
        {
            _id = glCreateProgram();
        }

        GLuint old_size = (GLuint)_shader.size();

        _shader.push_back(Shader());

        if (!_shader[old_size].createLoadAndCompileFromSourceFile(
                type, file_name, logging_is_enabled, output))
        {
            _shader.resize(old_size);
            return GL_FALSE;
        }

        glAttachShader(_id, _shader[old_size]._id);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::attachExistingShader(GLuint index)
    {
        if (!_id || index >= _shader.size())
        {
            return GL_FALSE;
        }

        GLboolean already_attached = GL_FALSE;

        GLint attached_shader_count = 0;
        glGetProgramiv(_id, GL_ATTACHED_SHADERS, &attached_shader_count);

        if (attached_shader_count)
        {
            GLuint *attached_shaders = new GLuint[attached_shader_count];
            glGetAttachedShaders(_id, attached_shader_count, nullptr, attached_shaders);

            for (GLuint i = 0; i < (GLuint)attached_shader_count && !already_attached;
                 i++)
            {
                already_attached = (_shader[index]._id == attached_shaders[i]);
            }

            delete [] attached_shaders;
        }

        if (already_attached)
        {
            return GL_FALSE;
        }

        glAttachShader(_id, _shader[index]._id);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::detachExistingShader(GLuint index)
    {
        if (!_id || index >= _shader.size())
        {
            return GL_FALSE;
        }

        GLboolean attached = GL_FALSE;

        GLint attached_shader_count = 0;
        glGetProgramiv(_id, GL_ATTACHED_SHADERS, &attached_shader_count);

        if (attached_shader_count)
        {
            GLuint *attached_shaders = new GLuint[attached_shader_count];
            glGetAttachedShaders(_id, attached_shader_count, nullptr, attached_shaders);

            for (GLuint i = 0; i < (GLuint)attached_shader_count && !attached; i++)
            {
                attached = (_shader[index]._id == attached_shaders[i]);
            }

            delete [] attached_shaders;
        }

        if (!attached)
        {
            return GL_FALSE;
        }

        glDetachShader(_id, _shader[index]._id);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::linkAttachedShaders(
            GLboolean logging_is_enabled, ostream &output)
    {
        if (!_id)
        {
            if (logging_is_enabled)
            {
                output << "The shader program object does not exist!" << endl;
            }

            return GL_FALSE;
        }

        if (logging_is_enabled)
        {
            output << "Linking the attached shaders (" << _shader.size() << ")...";
        }

        glLinkProgram(_id);
        glGetProgramiv(_id, GL_LINK_STATUS, &_link_status);

        if (!_link_status)
        {
            if (logging_is_enabled)
            {
                GLint info_log_length = 0;
                glGetProgramiv(_id, GL_INFO_LOG_LENGTH, &info_log_length);

                if (info_log_length > 0)
                {
                    GLchar *info_log = new GLchar[info_log_length];
                    glGetProgramInfoLog(_id, info_log_length, nullptr, info_log);

                    output << "\t\\begin{program information log}" << endl
                           << "\t\tid = " << _id << endl;
                    output <<  "\t\t" << info_log << endl;
                    output << "\t\\end{program information log}" << endl << endl;

                    delete[] info_log;
                }

                output << "\tUnsuccessful." << endl << endl;
            }

            return GL_FALSE;
        }

        if (logging_is_enabled)
        {
            output << "\tsuccessful." << endl << endl;
        }

        // loading the names and locations of active uniform variables
        if (logging_is_enabled)
        {
            output << "List of active uniform variables:" << endl;
        }

        _uniform.clear();

        GLint uniform_count;
        glGetProgramiv(_id, GL_ACTIVE_UNIFORMS, &uniform_count);

        GLsizei max_uniform_name_length;
        glGetProgramiv(_id, GL_ACTIVE_UNIFORM_MAX_LENGTH, &max_uniform_name_length);

        GLchar *uniform_name_data = new GLchar[max_uniform_name_length];

        for (GLint uniform = 0; uniform < uniform_count; uniform++)
        {
            GLsizei actual_length = 0;
            GLint   array_size    = 0;
            GLenum  type          = 0;

            glGetActiveUniform(_id, uniform, max_uniform_name_length, &actual_length,
                               &array_size, &type, uniform_name_data);
            string name(&uniform_name_data[0], actual_length);
            GLint  location = glGetUniformLocation(_id, name.c_str());

            if (logging_is_enabled)
            {
                output << "\tuniform: " << name
                       << ", type: " << openGLTypeToString(type)
                       << ", array size: " << array_size
                       << ", location: " << location << endl;
            }

            if (location >= 0)
            {
                _uniform[name] = Uniform(type, array_size, location);
            }
        }

        delete[] uniform_name_data;

        // loading the names and locations of active attributes
        if (logging_is_enabled)
        {
            output << "List of active attributes:" << endl;
        }

        _attribute.clear();

        GLint attribute_count;
        glGetProgramiv(_id, GL_ACTIVE_ATTRIBUTES, &attribute_count);

        GLsizei max_attribute_name_length;
        glGetProgramiv(_id, GL_ACTIVE_ATTRIBUTE_MAX_LENGTH, &max_attribute_name_length);

        GLchar *attribute_name_data = new GLchar[max_attribute_name_length];

        for (GLint attribute = 0; attribute < attribute_count; attribute++)
        {
            GLsizei actual_length = 0;
            GLint   array_size    = 0;
            GLenum  type          = 0;

            glGetActiveAttrib(_id, attribute, max_attribute_name_length, &actual_length,
                              &array_size, &type, attribute_name_data);
            string name(&attribute_name_data[0], actual_length);
            GLint  location = glGetAttribLocation(_id, name.c_str());

            if (logging_is_enabled)
            {
                output << "\tattribute: " << name
                       << ", type: " << openGLTypeToString(type)
                       << ", array size: " << array_size
                       << ", location: " << location << endl;
            }

            if (location >= 0)
            {
                _attribute[name] = Attribute(type, array_size, location);
            }
        }

        delete[] attribute_name_data;

        return GL_TRUE;
    }

    GLboolean ShaderProgram::enable(
            GLboolean logging_is_enabled, std::ostream &output) const
    {
        GLboolean all_shaders_were_compiled_successfully = GL_TRUE;

        for (GLuint i = 0; i < _shader.size() && all_shaders_were_compiled_successfully;
             i++)
        {
            all_shaders_were_compiled_successfully &= _shader[i]._compile_status;
        }

        if (_id && _link_status && all_shaders_were_compiled_successfully)
        {
            glUseProgram(_id);

            if (logging_is_enabled)
            {
                output << "Validating the program...";
            }

            glValidateProgram(_id);

            GLint validate_status = 0;
            glGetProgramiv(_id, GL_VALIDATE_STATUS, &validate_status);

            if (!validate_status)
            {
                if (logging_is_enabled)
                {
                    GLint info_log_length = 0;
                    glGetProgramiv(_id, GL_INFO_LOG_LENGTH, &info_log_length);

                    if (info_log_length > 0)
                    {
                        GLchar *info_log = new GLchar[info_log_length];
                        glGetProgramInfoLog(_id, info_log_length, nullptr, info_log);

                        output << "\t\\begin{program validation information log}" << endl
                               << "\t\tid = " << _id << endl;
                        output <<  "\t\t" << info_log << endl;
                        output << "\t\\end{program validation information log}" << endl
                               << endl;

                        delete[] info_log;

                        output << "\tUnsuccessful." << endl << endl;
                    }
                }

                return GL_FALSE;
            }

            if (logging_is_enabled)
            {
                output << "\tsuccessful." << endl << endl;
            }

            return GL_TRUE;
        }

        return GL_FALSE;
    }

    GLvoid ShaderProgram::disable() const
    {
        glUseProgram(0);
    }

    //(*@\Green{// getters}@*)

    GLint ShaderProgram::positionAttributeLocation() const
    {
        map<string, Attribute>::const_iterator iterator =
                _attribute.find(_position_attribute_name);

        if (iterator == _attribute.end() || iterator->second.type() != GL_FLOAT_VEC3)
        {
            return -1;
        }

        return iterator->second.location();
    }

    GLint ShaderProgram::normalAttributeLocation() const
    {
        map<string, Attribute>::const_iterator iterator =
                _attribute.find(_normal_attribute_name);

        if (iterator == _attribute.end() || iterator->second.type() != GL_FLOAT_VEC3)
        {
            return -1;
        }

        return iterator->second.location();
    }

    GLint ShaderProgram::colorAttributeLocation() const
    {
        map<string, Attribute>::const_iterator iterator =
                _attribute.find(_color_attribute_name);

        if (iterator == _attribute.end() || iterator->second.type() != GL_FLOAT_VEC4)
        {
            return -1;
        }

        return iterator->second.location();
    }

    GLint ShaderProgram::textureAttributeLocation() const
    {
        map<string, Attribute>::const_iterator iterator =
                _attribute.find(_texture_attribute_name);

        if (iterator == _attribute.end() || iterator->second.type() != GL_FLOAT_VEC4)
        {
            return -1;
        }

        return iterator->second.location();
    }

    ShaderProgram::Attribute ShaderProgram::activeAttribute(const string &name) const
    {
        map<string, Attribute>::const_iterator iterator = _attribute.find(name);

        if (iterator == _attribute.end())
        {
            return Attribute();
        }

        return iterator->second;
    }

    ShaderProgram::Uniform ShaderProgram::activeUniform(const string &name) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return Uniform();
        }

        return iterator->second;
    }

    //(*@\Green{// setters}@*)

    GLvoid ShaderProgram::setPositionAttributeName(const string &name)
    {
        _position_attribute_name = name;
    }

    GLvoid ShaderProgram::setNormalAttributeName(const string &name)
    {
        _normal_attribute_name = name;
    }

    GLvoid ShaderProgram::setColorAttributeName(const string &name)
    {
        _color_attribute_name = name;
    }

    GLvoid ShaderProgram::setTextureAttributeName(const string &name)
    {
        _texture_attribute_name = name;
    }

    GLboolean ShaderProgram::setUniformValue1i(const string &name, GLint value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform1i(iterator->second.location(), value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue2i(
        const string &name, GLint value_0, GLint value_1) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform2i(iterator->second.location(), value_0, value_1);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue3i(
        const string &name, GLint value_0, GLint value_1, GLint value2) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform3i(iterator->second.location(), value_0, value_1, value2);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue4i(
        const string &name,
        GLint value_0, GLint value_1, GLint value2, GLint value_3) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform4i(iterator->second.location(), value_0, value_1, value2, value_3);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue1ui(const string &name, GLuint value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform1ui(iterator->second.location(), value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue2ui(
        const string &name, GLuint value_0, GLuint value_1) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform2ui(iterator->second.location(), value_0, value_1);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue3ui(
        const string &name, GLuint value_0, GLuint value_1, GLuint value2) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform3ui(iterator->second.location(), value_0, value_1, value2);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue4ui(
        const string &name,
        GLuint value_0, GLuint value_1, GLuint value2, GLuint value_3) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform4ui(iterator->second.location(), value_0, value_1, value2, value_3);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue1f(const string &name, GLfloat value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform1f(iterator->second.location(), value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue2f(
        const string &name, GLfloat value_0, GLfloat value_1) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform2f(iterator->second.location(), value_0, value_1);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue3f(
        const string &name, GLfloat value_0, GLfloat value_1, GLfloat value2) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform3f(iterator->second.location(), value_0, value_1, value2);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue4f(
        const string &name,
        GLfloat value_0, GLfloat value_1, GLfloat value2, GLfloat value_3) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform4f(iterator->second.location(), value_0, value_1, value2, value_3);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue1iv(
        const string &name, GLsizei count, const GLint *value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform1iv(iterator->second.location(), count, value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue2iv(
        const string &name, GLsizei count, const GLint *value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform2iv(iterator->second.location(), count, value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue3iv(
        const string &name, GLsizei count, const GLint *value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform3iv(iterator->second.location(), count, value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue4iv(
        const string &name, GLsizei count, const GLint *value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform4iv(iterator->second.location(), count, value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue1uiv(
        const string &name, GLsizei count, const GLuint *value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform1uiv(iterator->second.location(), count, value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue2uiv(
        const string &name, GLsizei count, const GLuint *value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform2uiv(iterator->second.location(), count, value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue3uiv(
        const string &name, GLsizei count, const GLuint *value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform3uiv(iterator->second.location(), count, value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue4uiv(
        const string &name, GLsizei count, const GLuint *value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform4uiv(iterator->second.location(), count, value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue1fv(
        const string &name, GLsizei count, const GLfloat *value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform1fv(iterator->second.location(), count, value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue2fv(
        const string &name, GLsizei count, const GLfloat *value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform2fv(iterator->second.location(), count, value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue3fv(
        const string &name, GLsizei count, const GLfloat *value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform3fv(iterator->second.location(), count, value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformValue4fv(
        const string &name, GLsizei count, const GLfloat *value) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform4fv(iterator->second.location(), count, value);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformMatrix2fv(
        const std::string &name, GLsizei count,
        GLboolean transpose, const GLfloat *values) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniformMatrix2fv(iterator->second.location(), count, transpose, values);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformMatrix3fv(
        const std::string &name, GLsizei count,
        GLboolean transpose, const GLfloat *values) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniformMatrix3fv(iterator->second.location(), count, transpose, values);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformMatrix4fv(
        const std::string &name, GLsizei count,
        GLboolean transpose, const GLfloat *values) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniformMatrix4fv(iterator->second.location(), count, transpose, values);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformMatrix2x3fv(
        const std::string &name, GLsizei count,
        GLboolean transpose, const GLfloat *values) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniformMatrix2x3fv(iterator->second.location(), count, transpose, values);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformMatrix3x2fv(
        const std::string &name, GLsizei count,
        GLboolean transpose, const GLfloat *values) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniformMatrix3x2fv(iterator->second.location(), count, transpose, values);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformMatrix2x4fv(
        const std::string &name, GLsizei count,
        GLboolean transpose, const GLfloat *values) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniformMatrix2x4fv(iterator->second.location(), count, transpose, values);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformMatrix4x2fv(
        const std::string &name, GLsizei count,
        GLboolean transpose, const GLfloat *values) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniformMatrix4x2fv(iterator->second.location(), count, transpose, values);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformMatrix3x4fv(
        const std::string &name, GLsizei count,
        GLboolean transpose, const GLfloat *values) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniformMatrix3x4fv(iterator->second.location(), count, transpose, values);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformMatrix4x3fv(
        const std::string &name, GLsizei count,
        GLboolean transpose, const GLfloat *values) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniformMatrix4x3fv(iterator->second.location(), count, transpose, values);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformColor(
        const std::string &name, const Color4 &color) const
    {
        map<string, Uniform>::const_iterator iterator = _uniform.find(name);

        if (iterator == _uniform.end())
        {
            return GL_FALSE;
        }

        glUniform4fv(iterator->second.location(), 1, color.address());

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformDirectionalLight(
        const string &name, const DirectionalLight &light) const
    {
        map<string, Uniform>::const_iterator position_iterator    =
                _uniform.find(name + ".position");
        map<string, Uniform>::const_iterator half_vector_iterator =
                _uniform.find(name + ".half_vector");
        map<string, Uniform>::const_iterator ambient_iterator     =
                _uniform.find(name + ".ambient");
        map<string, Uniform>::const_iterator diffuse_iterator     =
                _uniform.find(name + ".diffuse");
        map<string, Uniform>::const_iterator specular_iterator    =
                _uniform.find(name + ".specular");

        map<string, Uniform>::const_iterator end = _uniform.end();

        if (position_iterator == end || half_vector_iterator == end ||
            ambient_iterator  == end || diffuse_iterator     == end ||
            specular_iterator == end)
        {
            return GL_FALSE;
        }

        glUniform4fv(position_iterator->second.location(), 1,
                     light.addressOfPosition());
        glUniform4fv(half_vector_iterator->second.location(), 1,
                     light.addressOfHalfVector());
        glUniform4fv(ambient_iterator->second.location(), 1,
                     light.addressOfAmbientIntensity());
        glUniform4fv(diffuse_iterator->second.location(), 1,
                     light.addressOfDiffuseIntensity());
        glUniform4fv(specular_iterator->second.location(), 1,
                     light.addressOfSpecularIntensity());

        map<string, Uniform>::const_iterator spot_cos_cutoff_iterator       =
                _uniform.find(name + ".spot_cos_cutoff");
        map<string, Uniform>::const_iterator constant_attenuation_iterator  =
                _uniform.find(name + ".constant_attenuation");
        map<string, Uniform>::const_iterator linear_attenuation_iterator    =
                _uniform.find(name + ".linear_attenuation");
        map<string, Uniform>::const_iterator quadratic_attenuation_iterator =
                _uniform.find(name + ".quadratic_attenuation");
        map<string, Uniform>::const_iterator spot_direction_iterator        =
                _uniform.find(name + ".spot_direction");
        map<string, Uniform>::const_iterator spot_exponent_iterator         =
                _uniform.find(name + ".spot_exponent");

        if (spot_cos_cutoff_iterator    == end || constant_attenuation_iterator  == end ||
            linear_attenuation_iterator == end || quadratic_attenuation_iterator == end ||
            spot_direction_iterator     == end || spot_exponent_iterator         == end)
        {
            return GL_FALSE;
        }

        glUniform1f(spot_cos_cutoff_iterator->second.location(), 180.0f);
        glUniform1f(constant_attenuation_iterator->second.location(), 1.0f);
        glUniform1f(linear_attenuation_iterator->second.location(), 0.0f);
        glUniform1f(quadratic_attenuation_iterator->second.location(), 0.0f);
        glUniform3f(spot_direction_iterator->second.location(), 0.0f, 0.0f, -1.0f);
        glUniform1f(spot_exponent_iterator->second.location(), 0.0);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformPointLight(
        const string &name, const PointLight &light) const
    {
        map<string, Uniform>::const_iterator position_iterator              =
                _uniform.find(name + ".position");
        map<string, Uniform>::const_iterator half_vector_iterator           =
                _uniform.find(name + ".half_vector");
        map<string, Uniform>::const_iterator ambient_iterator               =
                _uniform.find(name + ".ambient");
        map<string, Uniform>::const_iterator diffuse_iterator               =
                _uniform.find(name + ".diffuse");
        map<string, Uniform>::const_iterator specular_iterator              =
                _uniform.find(name + ".specular");
        map<string, Uniform>::const_iterator spot_cos_cutoff_iterator       =
                _uniform.find(name + ".spot_cos_cutoff");
        map<string, Uniform>::const_iterator constant_attenuation_iterator  =
                _uniform.find(name + ".constant_attenuation");
        map<string, Uniform>::const_iterator linear_attenuation_iterator    =
                _uniform.find(name + ".linear_attenuation");
        map<string, Uniform>::const_iterator quadratic_attenuation_iterator =
                _uniform.find(name + ".quadratic_attenuation");

        map<string, Uniform>::const_iterator end = _uniform.end();

        if (position_iterator           == end || half_vector_iterator           == end ||
            ambient_iterator            == end || diffuse_iterator               == end ||
            specular_iterator           == end ||
            spot_cos_cutoff_iterator    == end || constant_attenuation_iterator  == end ||
            linear_attenuation_iterator == end || quadratic_attenuation_iterator == end)
        {
            return GL_FALSE;
        }

        glUniform4fv(position_iterator->second.location(), 1,
                     light.addressOfPosition());
        glUniform4fv(half_vector_iterator->second.location(), 1,
                     light.addressOfHalfVector());
        glUniform4fv(ambient_iterator->second.location(), 1,
                     light.addressOfAmbientIntensity());
        glUniform4fv(diffuse_iterator->second.location(), 1,
                     light.addressOfDiffuseIntensity());
        glUniform4fv(specular_iterator->second.location(), 1,
                     light.addressOfSpecularIntensity());

        glUniform1f(spot_cos_cutoff_iterator->second.location(), 180.0f);
        glUniform1f(constant_attenuation_iterator->second.location(),
                    light.constantAttenuation());
        glUniform1f(linear_attenuation_iterator->second.location(),
                    light.linearAttenuation());
        glUniform1f(quadratic_attenuation_iterator->second.location(),
                    light.quadraticAttenuation());

        map<string, Uniform>::const_iterator spot_direction_iterator =
                _uniform.find(name + ".spot_direction");
        map<string, Uniform>::const_iterator spot_exponent_iterator  =
                _uniform.find(name + ".spot_exponent");

        if (spot_direction_iterator == end || spot_exponent_iterator == end)
        {
            return GL_FALSE;
        }

        glUniform3f(spot_direction_iterator->second.location(), 0.0f, 0.0f, -1.0f);
        glUniform1f( spot_exponent_iterator->second.location(), 0.0f);

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformSpotlight(
        const string &name, const Spotlight &light) const
    {
        map<string, Uniform>::const_iterator position_iterator              =
                _uniform.find(name + ".position");
        map<string, Uniform>::const_iterator half_vector_iterator           =
                _uniform.find(name + ".half_vector");
        map<string, Uniform>::const_iterator ambient_iterator               =
                _uniform.find(name + ".ambient");
        map<string, Uniform>::const_iterator diffuse_iterator               =
                _uniform.find(name + ".diffuse");
        map<string, Uniform>::const_iterator specular_iterator              =
                _uniform.find(name + ".specular");
        map<string, Uniform>::const_iterator spot_cos_cutoff_iterator       =
                _uniform.find(name + ".spot_cos_cutoff");
        map<string, Uniform>::const_iterator constant_attenuation_iterator  =
                _uniform.find(name + ".constant_attenuation");
        map<string, Uniform>::const_iterator linear_attenuation_iterator    =
                _uniform.find(name + ".linear_attenuation");
        map<string, Uniform>::const_iterator quadratic_attenuation_iterator =
                _uniform.find(name + ".quadratic_attenuation");
        map<string, Uniform>::const_iterator spot_direction_iterator        =
                _uniform.find(name + ".spot_direction");
        map<string, Uniform>::const_iterator spot_exponent_iterator         =
                _uniform.find(name + ".spot_exponent");

        map<string, Uniform>::const_iterator end = _uniform.end();

        if (position_iterator           == end || half_vector_iterator           == end ||
            ambient_iterator            == end || diffuse_iterator               == end ||
            specular_iterator           == end ||
            spot_cos_cutoff_iterator    == end || constant_attenuation_iterator  == end ||
            linear_attenuation_iterator == end || quadratic_attenuation_iterator == end ||
            spot_direction_iterator     == end || spot_exponent_iterator         == end)
        {
            return GL_FALSE;
        }

        const GLdouble *spot_direction = light.addressOfSpotDirection();

        glUniform4fv(position_iterator->second.location(), 1,
                     light.addressOfPosition());
        glUniform4fv(half_vector_iterator->second.location(), 1,
                     light.addressOfHalfVector());
        glUniform4fv(ambient_iterator->second.location(), 1,
                     light.addressOfAmbientIntensity());
        glUniform4fv(diffuse_iterator->second.location(), 1,
                     light.addressOfDiffuseIntensity());
        glUniform4fv(specular_iterator->second.location(), 1,
                     light.addressOfSpecularIntensity());
        glUniform1f (spot_cos_cutoff_iterator->second.location(),
                     light.spotCosCutoff());
        glUniform1f (constant_attenuation_iterator->second.location(),
                     light.constantAttenuation());
        glUniform1f (linear_attenuation_iterator->second.location(),
                     light.linearAttenuation());
        glUniform1f (quadratic_attenuation_iterator->second.location(),
                     light.quadraticAttenuation());
        glUniform3f (spot_direction_iterator->second.location(),
                    (GLfloat)spot_direction[0], (GLfloat)spot_direction[1],
                    (GLfloat)spot_direction[2]);
        glUniform1f (spot_exponent_iterator->second.location(), light.spotExponent());

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformMaterial(
        const std::string &name, const Material &material) const
    {
        map<string, Uniform>::const_iterator ambient_iterator   =
                _uniform.find(name + ".ambient");
        map<string, Uniform>::const_iterator diffuse_iterator   =
                _uniform.find(name + ".diffuse");
        map<string, Uniform>::const_iterator specular_iterator  =
                _uniform.find(name + ".specular");
        map<string, Uniform>::const_iterator emission_iterator  =
                _uniform.find(name + ".emission");
        map<string, Uniform>::const_iterator shininess_iterator =
                _uniform.find(name + ".shininess");

        map<string, Uniform>::const_iterator end = _uniform.end();

        if (ambient_iterator  == end || diffuse_iterator  == end ||
            specular_iterator == end || emission_iterator == end ||
            shininess_iterator == end)
        {
            return GL_FALSE;
        }

        glUniform4fv(ambient_iterator->second.location(), 1,
                     material.addressOfAmbientReflectionCoefficients());
        glUniform4fv(diffuse_iterator->second.location(), 1,
                     material.addressOfDiffuseReflectionCoefficients());
        glUniform4fv(specular_iterator->second.location(), 1,
                     material.addressOfSpecularReflectionCoefficients());
        glUniform4fv(emission_iterator->second.location(), 1,
                     material.addressOfEmissionColor());
        glUniform1f (shininess_iterator->second.location(),
                     material.shininess());

        return GL_TRUE;
    }

    GLboolean ShaderProgram::setUniformColorMaterial(
        const std::string &name,
        const Color4 &ambient_and_diffuse_reflection_coefficients,
        const Color4 &specular_reflection_coefficients,
        GLfloat shininess) const
    {
        map<string, Uniform>::const_iterator ambient_iterator   =
                _uniform.find(name + ".ambient");
        map<string, Uniform>::const_iterator diffuse_iterator   =
                _uniform.find(name + ".diffuse");
        map<string, Uniform>::const_iterator specular_iterator  =
                _uniform.find(name + ".specular");
        map<string, Uniform>::const_iterator emission_iterator  =
                _uniform.find(name + ".emission");
        map<string, Uniform>::const_iterator shininess_iterator =
                _uniform.find(name + ".shininess");

        map<string, Uniform>::const_iterator end = _uniform.end();

        if (ambient_iterator  == end || diffuse_iterator  == end ||
            specular_iterator == end || emission_iterator == end ||
            shininess_iterator == end)
        {
            return GL_FALSE;
        }

        glUniform4fv(ambient_iterator->second.location(), 1,
                     ambient_and_diffuse_reflection_coefficients.address());
        glUniform4fv(diffuse_iterator->second.location(), 1,
                     ambient_and_diffuse_reflection_coefficients.address());
        glUniform4fv(specular_iterator->second.location(), 1,
                     specular_reflection_coefficients.address());
        glUniform4f (emission_iterator->second.location(), 0.0f, 0.0f, 0.0f, 0.0f);
        glUniform1f (shininess_iterator->second.location(), shininess);

        return GL_TRUE;
    }

    //(*@\Green{// destructor}@*)
    ShaderProgram::~ShaderProgram()
    {
        for (GLuint i = 0; i < _shader.size(); i++)
        {
            glDetachShader(_id, _shader[i]._id);
            _shader[i].flagForDeletion();
        }

        if (_id)
        {
            glDeleteProgram(_id);
        }
    }
}
