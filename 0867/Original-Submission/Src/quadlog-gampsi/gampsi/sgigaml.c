/* Fortran dgamma() using SGI IRIX C extended precision library */

#include <math.h>

double dgamma_(double *x)
{
    long double lg;
    long double xx;
    xx = *x;
    lg = lgammal(xx);

    return ((double)signgaml*expl(lg));
}
