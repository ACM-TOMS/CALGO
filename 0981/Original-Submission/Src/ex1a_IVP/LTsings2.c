/**     EXAMPLE 1a

    Application of Talbot's method to solve
    the following PDE problem

            u_t (x,t) = u_x (x,t),  x>X0,   t>0
            u(x,0+) = x
            u(X0,t) = X0 + t

    whose analytical solution is

            u(x,t) = x + t

    and Laplace Transform is

            U(x,s) = x/s + 1/s^2

    where s = 0 is a double pole.
 */

#include <complex.h>


unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0)
/* Abscissa of convergence and singularities of F(s) with their polar multiplicities */
{
    unsigned int NsingsTOT = 1; /* total number of singularities s_j */
    *Nsings = 1;                /* number of singularities s_j with Im(s_j) ge 0 */
    *sigma0 = 0;                /* abscissa of convergence */

    SINGS[0]=0.0+0.0*I;
    MULT[0]=2;

    return NsingsTOT;
}
