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
