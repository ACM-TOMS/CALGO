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

void main(void)
{
    interpolated_normal   = normalize(vec3(N * vec4(normal, 0.0)));
    interpolated_position = vec3(VM * vec4(position, 1.0));
    interpolated_color    = color;
    interpolated_texture  = texture;
    gl_Position           = PVM * vec4(position, 1.0);
}
