//----------------------------------------------------------------------------------
// File:        Shaders/two_sided_lighting.frag
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

in vec3 interpolated_position;
in vec3 interpolated_normal;
in vec4 interpolated_color;
in vec4 interpolated_texture;

out vec4 fragment_color;

uniform float     transparency = 0.0;
uniform sampler2D sampler_2D_texture;

struct Material
{
    vec4    ambient;
    vec4    diffuse;
    vec4    specular;
    vec4    emission;
    float   shininess;
};

uniform Material front_material;
uniform Material back_material;

struct LightSource
{
    bool    enabled;
    vec4    position;
    vec4    half_vector;
    vec4    ambient;
    vec4    diffuse;
    vec4    specular;
    float   spot_cos_cutoff;
    float   constant_attenuation;
    float   linear_attenuation;
    float   quadratic_attenuation;
    vec3    spot_direction;
    float   spot_exponent;
};

const int number_of_light_sources = 1;

uniform LightSource light_source[number_of_light_sources];

// distance is measured from current vertex to the ith light source.
float calculateAttenuation(in int i, in float distance)
{
    return(1.0 / (light_source[i].constant_attenuation +
                  light_source[i].linear_attenuation * distance +
                  light_source[i].quadratic_attenuation * distance * distance));
}

// N denotes the unit varying normal vector.
void directionalLight(in int i, in vec3 N, in float shininess,
                      inout vec4 ambient, inout vec4 diffuse, inout vec4 specular)
{
    vec3 L = normalize(light_source[i].position.xyz);
   
    float N_dot_L = dot(N, L);
   
    if (N_dot_L > 0.0)
    {   
        vec3 H = light_source[i].half_vector.xyz;
       
        float pf = pow(max(dot(N, H), 0.0), shininess);

        diffuse  += light_source[i].diffuse  * N_dot_L;
        specular += light_source[i].specular * pf;
        ambient  += light_source[i].ambient;
    }
   
}

// N denotes the unit varying normal vector, while V corresponds to the varying vertex position.
void pointLight(in int i, in vec3 N, in vec3 V, in float shininess,
                inout vec4 ambient, inout vec4 diffuse, inout vec4 specular)
{
    vec3 D = light_source[i].position.xyz - V;
    vec3 L = normalize(D);

    float distance = length(D);
    float attenuation = calculateAttenuation(i, distance);

    float N_dot_L = dot(N, L);

    if (N_dot_L > 0.0)
    {   
        vec3 E = normalize(-V);
        vec3 R = reflect(-L, N);
       
        float pf = pow(max(dot(R,E), 0.0), shininess);

        diffuse  += light_source[i].diffuse  * attenuation * N_dot_L;
        specular += light_source[i].specular * attenuation * pf;
        ambient  += light_source[i].ambient  * attenuation;
    }
   
}

// N denotes the unit varying normal vector, while V corresponds to the varying vertex position.
void spotlight(in int i, in vec3 N, in vec3 V, in float shininess,
               inout vec4 ambient, inout vec4 diffuse, inout vec4 specular)
{
    vec3 D = light_source[i].position.xyz - V;
    vec3 L = normalize(D);

    float distance = length(D);
    float attenuation = calculateAttenuation(i, distance);

    float N_dot_L = dot(N,L);

    if (N_dot_L > 0.0)
    {   
        float spot_effect = dot(normalize(light_source[i].spot_direction), -L);
       
        if (spot_effect > light_source[i].spot_cos_cutoff)
        {
            attenuation *=  pow(spot_effect, light_source[i].spot_exponent);

            vec3 E = normalize(-V);
            vec3 R = reflect(-L, N);
       
            float pf = pow(max(dot(R, E), 0.0), shininess);

            diffuse  += light_source[i].diffuse  * attenuation * N_dot_L;
            specular += light_source[i].specular * attenuation * pf;
        }
        ambient  += light_source[i].ambient * attenuation;
    }
   
}

// N denotes the unit varying normal vector, while V corresponds to the varying vertex position.
void calculateLighting(in int number_of_light_sources, in vec3 N, in vec3 V, in float shininess,
                       inout vec4 ambient, inout vec4 diffuse, inout vec4 specular)
{
    // Just loop through each light, and if its enabled add
    // its contributions to the color of the pixel.
    for (int i = 0; i < number_of_light_sources; i++)
    {
        if (light_source[i].enabled)
        {
            if (light_source[i].position.w == 0.0)
            {
                directionalLight(i, N, shininess, ambient, diffuse, specular);
            }
            else
            {
                if (light_source[i].spot_cos_cutoff == 180.0)
                {
                    pointLight(i, N, V, shininess, ambient, diffuse, specular);
                }
                else
                {
                    spotlight(i, N, V, shininess, ambient, diffuse, specular);
                }
            }
        }
    }
}

void main()
{
    // Normalize the interpolated normal. A varying variable CANNOT be modified by a fragment shader.
    // So a new variable needs to be created.
    vec3 n = normalize(interpolated_normal);

    vec4 ambient = vec4(0.0), diffuse = vec4(0.0), specular = vec4(0.0), color = vec4(0.0);

    // Initialize the contributions for the front-face-pass over the lights.
    ambient  = texture2DProj(sampler_2D_texture, interpolated_texture);
    diffuse  = texture2DProj(sampler_2D_texture, interpolated_texture);

    calculateLighting(number_of_light_sources, n, interpolated_position, front_material.shininess,
                      ambient, diffuse, specular);

    color += front_material.emission +
             (ambient  * (front_material.ambient + interpolated_color)) +
             (diffuse  * (front_material.diffuse + interpolated_color)) +
             (specular * (front_material.specular));

    // Re-initialize the contributions for the back-face-pass over the lights.
    ambient  = texture2DProj(sampler_2D_texture, interpolated_texture);
    diffuse  = texture2DProj(sampler_2D_texture, interpolated_texture);
    specular = vec4(0.0);

    // Now caculate the back contribution. All that needs to be done is to flip the normal.
    calculateLighting(number_of_light_sources, -n, interpolated_position, back_material.shininess,
                      ambient, diffuse, specular);

    color += back_material.emission +
             (ambient  * (back_material.ambient + interpolated_color)) +
             (diffuse  * (back_material.diffuse + interpolated_color)) +
             (specular * (back_material.specular));

    color = clamp(color, 0.0, 1.0);

    fragment_color = vec4(color.rgb, clamp(1.0 - transparency, 0.0, 1.0));
}
