/* Fortran dgamma() using GNU Scientic Library (gsl) */

#include <gsl_errno.h>
#include <gsl_sf_gamma.h>

double dgamma_(double *x)
{
    gsl_sf_result result;
    if (gsl_sf_gamma_e(*x, &result) == GSL_SUCCESS)
	return (result.val);
    else
    {
	result.val = 0.0;
	return (result.val / result.val);
    }
}

double dgamma(double *x)
{
    return (dgamma_(x));    
}
