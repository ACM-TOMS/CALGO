#version 130

uniform vec4 color;          (*@\Green{// uniform color variable that can be modified by the user\label{src:color.frag:color}}@*)
out     vec4 fragment_color; (*@\Green{// output of the current fragment shader}@*)

void main()
{
    fragment_color = color;
}
