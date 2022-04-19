/**      EXAMPLE 0

    Laplace Transform function (LTfun) [F24(s) from test functions]:
        F(s) = F(s) = s/(s^2+9)^2

    2 polar singularities (double poles):
        s_j = {+/-3*i}

    abscissa of convergence:
        sigma0 = 0

    Inverse Laplace Transform function (ILTfun):
        f(t) = t*sin(3*t)/6
 */

#include <complex.h>


unsigned int LTsings(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0)
/* Abscissa of convergence and singularities of F(s) with their polar multiplicities */
{   unsigned int NsingsTOT = 2; /* total number of singularities s_j */
    *Nsings = 1;                /* number of singularities s_j with Im(s_j) ge 0 */
    *sigma0 = 0;                /* abscissa of convergence */

    SINGS[0]=0.0+3.0*I; SINGS[1]=0.0-3.0*I;
    MULT[0]=2; MULT[1]=2;

    return NsingsTOT;
}
