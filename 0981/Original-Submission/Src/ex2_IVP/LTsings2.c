/**	        Example 2

    Application of Talbot's method to solve
    the following PDE problem

            u_tx (x,t) = exp(-x)*cos(t),    x>0,    t>0
            u_x(x,0+)  = 0
            u(0,t)     = 0

    The analytical solution is

            u(x,t) = sin(t)*(1-exp(-x))

    and Laplace Transform is

            U(x,s) = (1-exp(-x))/(s^2+1)

    where s = +/-i are simple poles.

    The ODE problem is

       U' = exp(-x)/(s^2+1),        x>0
       U(0) = 0
*/

#include <complex.h>


unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0)
/* Abscissa of convergence and singularities of F(s) with their polar multiplicities */
{   unsigned int NsingsTOT = 2; /* total number of singularities s_j */
    *Nsings = 1;                /* number of singularities s_j with Im(s_j) ge 0 */
    *sigma0 = 0;                /* abscissa of convergence */

    SINGS[0]=0.0 + 1.0*I;
    MULT[0]=1;
    SINGS[1]=0.0 - 1.0*I;
    MULT[1]=1;

    return NsingsTOT;
}
