/**     EXAMPLE 3a

    Application of Talbot's method to solve
    the following PDE problem

        u_t (x,t) = u_xx (x,t),		x>0, t>0
		  u(x,0+) = x*(x-1)
		  u(0,t)  = 2*t
		  u_x(0,t)= -1

	in [0,X1]x[T0,T1].

    The analytical solution is

            u(x,t) = x*(x-1) + 2*t

    and Laplace Transform is

                      2    x*(x-1)
            U(x,s) = --- + -------
                     s^2      s

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
