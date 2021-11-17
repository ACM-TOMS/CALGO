#version 130

uniform mat4 PVM;       (*@\Green{// product of projection, view and model matrices\label{src:color.vert:PVM}}@*)
in      vec3 position;  (*@\Green{// attribute}@*)

void main()
{
    (*@\Green{// convert the given position to clip coordinates and pass along}@*)
    gl_Position = PVM * vec4(position, 1.0);
}
