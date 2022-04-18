(*@\Green{// \ldots}@*)
(*@\Green{// these lines coincide with lines \mref{src:two_sided_lighting.frag:part_1:start}--\mref{src:two_sided_lighting.frag:part_1:end} of Listing \mref{src:two_sided_lighting.frag}}@*)
(*@\Green{// \ldots}\Suppressnumber{}@*)
(*@\Red{\rule{\textwidth}{0.75pt}}\Reactivatenumber{}@*)
uniform float     scale_factor;          (*@\Green{// affects the number of reflection lines, its value should be}\Suppressnumber{}@*)
                                         (*@\Green{// greater than zero}\Reactivatenumber{}@*)
uniform float     smoothing;             (*@\Green{// smooths the common boundaries of adjacent stripes,}\Suppressnumber{}@*)
                                         (*@\Green{// its value should be greater than zero}\Reactivatenumber{}@*)
uniform float     shading;               (*@\Green{// shades the model, its value should be in the range $\left[0,1\right]$}@*)

const   float     TWO_PI = 6.283185307179586476925286766559;    (*@\Suppressnumber{}@*)
(*@\Blue{\rule{\textwidth}{0.75pt}}\Reactivatenumber{}@*)
(*@\Green{// \ldots}@*)
(*@\Green{// these lines coincide with lines \mref{src:two_sided_lighting.frag:part_2:start}--\mref{src:two_sided_lighting.frag:part_2:end} of Listing \mref{src:two_sided_lighting.frag}}@*)
(*@\Green{// \ldots}@*)
void main()
{
    (*@\Green{// \ldots}@*)
    (*@\Green{// these lines coincide with lines \mref{src:two_sided_lighting.frag:part_3:start}--\mref{src:two_sided_lighting.frag:part_3:end} of Listing \mref{src:two_sided_lighting.frag}}@*)
    (*@\Green{// \ldots}\Suppressnumber{}@*)
(*@\Red{\rule{\textwidth}{0.75pt}}\Reactivatenumber{}@*)
    (*@\Green{// generating reflection lines\ldots}@*)
    vec3 v  = normalize(interpolated_position);
    vec3 r  = v - 2.0 * n * dot(n, v); (*@\Green{// reflection vector}@*)
    r.z    += 1.0;
    r.z     = 0.5 / sqrt(dot(r, r));
    r.xy    = (r.xy * r.z) + 0.5;
    r      *= 2.0;

    float sharpness = 1.0 / asin(smoothing * fwidth(r.x) * 2.0 * scale_factor);

    color   *= vec4(clamp(0.5 + sharpness * sin(TWO_PI * r.x * scale_factor), 0.0, 1.0));
    color.a  = 1.0;

    fragment_color = clamp(dot(light_source[0].position.xyz, n) * shading +
                           (1.0 - shading), 0.0, 1.0) *
                     min(vec4(color.rgb, clamp(1.0 - transparency, 0.0, 1.0)), vec4(1.0)); (*@\Suppressnumber{}@*)
(*@\Blue{\rule{\textwidth}{0.75pt}}\Reactivatenumber{}@*)
}
