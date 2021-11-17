/* Fortran dgamma() using SGI IRIX C library */

#include <math.h>

double dgamma_(double *x)
{
    double lg;
    
    lg = lgamma(*x);
    return (signgam*exp(lg));
}

