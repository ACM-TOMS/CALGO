/* Fortran dgamma() using IBM RS/6000 AIX C library */

#include <math.h>

double dgamma(double *x)
{
    double lg;
    
    lg = lgamma(*x);
    return (signgam*exp(lg));
}

