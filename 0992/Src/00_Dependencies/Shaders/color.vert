//----------------------------------------------------------------------------------
// File:        Shaders/color.vert
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

uniform mat4 PVM;       // product of projection, view and model matrices
in      vec3 position;  // attribute

void main()
{
    // convert the given vertex to clip coordinates and pass along
    gl_Position = PVM * vec4(position, 1.0);
}
