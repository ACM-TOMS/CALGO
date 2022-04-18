#version 130

uniform mat4 VM;    (*@\Green{// product of view and model matrices}@*)
uniform mat4 PVM;   (*@\Green{// product of projection, view and model matrices}@*)
uniform mat4 N;     (*@\Green{// normal matrix, i.e., transposed inverse of VM}@*)

in vec3 position;   (*@\Green{// attributes associated with vertices}@*)
in vec3 normal;
in vec4 color;
in vec4 texture;

out vec3 interpolated_position; (*@\Green{// outputs of the current vertex shader (i.e., inputs of the fragment shader}@*)
out vec3 interpolated_normal;   (*@\Green{// listed in Listing \mref{src:two_sided_lighting.frag})}@*)
out vec4 interpolated_color;
out vec4 interpolated_texture;

void main()
{
    (*@\Green{// transform the normal vector to the eye space, then normalize it}@*)
    interpolated_normal   = normalize(vec3(N * vec4(normal, 0.0)));
   
    (*@\Green{// transform the vertex position to the eye space}@*)
    interpolated_position = vec3(VM * vec4(position, 1.0));

    interpolated_color    = color;
    interpolated_texture  = texture;
       
    (*@\Green{// convert the given position to clip coordinates and pass along}@*)
    gl_Position           = PVM * vec4(position, 1.0);
}
