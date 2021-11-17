//----------------------------------------------------------------------------------
// File:        Shaders/two_sided_lighting.vert
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 - August 15, 2017 - July 09, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
// Reference:   Randi J. Rost and Bill Licea-Kane. 2006. OpenGL shading language,
//              2nd edition. Addison-Wesley Professional, USA.
//----------------------------------------------------------------------------------

#version 130

uniform mat4 VM;    // product of view and model matrices
uniform mat4 PVM;   // product of projection, view and model matrices
uniform mat4 N;     // normal matrix, i.e., transposed inverse of VM

in vec3 position;
in vec3 normal;
in vec4 color;
in vec4 texture;

out vec3 interpolated_position;
out vec3 interpolated_normal;
out vec4 interpolated_color;
out vec4 interpolated_texture;

void main()
{
    // transform the normal vector to the eye space, then normalize it
    interpolated_normal   = normalize(vec3(N * vec4(normal, 0.0)));
   
    // transform the vertex position to the eye space
    interpolated_position = vec3(VM * vec4(position, 1.0));

    interpolated_color    = color;
    interpolated_texture  = texture;
       
    // convert the given vertex to clip coordinates and pass along
    gl_Position           = PVM * vec4(position, 1.0);
}
