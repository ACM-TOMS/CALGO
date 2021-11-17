/* Fortran dgamma() using Sun Solaris 2.x C library */

#include <math.h>

double dgamma_(double *x)
{
    double lg;
    
    lg = lgamma(*x);
    return (signgam*exp(lg));
}

