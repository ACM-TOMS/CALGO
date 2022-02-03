/* Fortran dgamma() using GNU gcc runtime library */

#include <math.h>

double dgamma_(double *x)
{
    long double lg;
    long double xx;
    xx = *x;
    lg = lgammal(xx);

    return ((double)signgam*expl(lg));
}
