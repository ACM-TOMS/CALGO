//----------------------------------------------------------------------------------
// File:        Shaders/color.frag
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

uniform vec4 color;
out     vec4 fragment_color;

void main()
{
    fragment_color = color;
}
