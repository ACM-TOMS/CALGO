/* Fortran dgamma() using IBM RS/6000 AIX C extended precision library */

#include <math.h>

double dgamma(double *x)
{
    long double lg;
    long double xx;
    xx = *x;
    lg = lgammal(xx);

    return ((double)signgam*expl(lg));
}
