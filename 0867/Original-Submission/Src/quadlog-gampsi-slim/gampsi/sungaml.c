/* Fortran dgamma() using Sun Solaris 2.x C extended precision library */

#include <sunmath.h>

double dgamma_(double *x)
{
    long double lg;
    long double xx;
    xx = *x;
    lg = lgammal(xx);

    return ((double)signgaml*expl(lg));
}
