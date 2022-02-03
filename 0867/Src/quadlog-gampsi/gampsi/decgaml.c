/* Fortran dgamma() using DEC Alpha OSF/1 C extended precision library */

#include <math.h>

double dgamma_(double *x)
{
    long double lg;
    long double xx;
    xx = *x;
    lg = lgammal(xx);

    return ((double)signgam*expl(lg));
}
