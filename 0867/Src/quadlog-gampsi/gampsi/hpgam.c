/* Fortran dgamma() using Hewlett-Packard HP-UX library */

#include <math.h>

double dgamma(double *x)
{
    double lg;
    
    lg = lgamma(*x);
    return (signgam*exp(lg));
}

