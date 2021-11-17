/* Fortran dpsi() using GNU Scientic Library (gsl) */

#include <gsl_errno.h>
#include <gsl_sf_psi.h>

double dpsi_(double *x)
{
    gsl_sf_result result;
    if (gsl_sf_psi_e(*x, &result) == GSL_SUCCESS)
	return (result.val);
    else
    {
	result.val = 0.0;
	return (result.val / result.val);
    }
}

double dpsi(double *x)
{
    return (dpsi_(x));    
}
