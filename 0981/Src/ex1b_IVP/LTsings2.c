/**     EXAMPLE 1b

    Application of Talbot's method to solve
    the following PDE problem

            u_t (x,t) = u_x (x,t),  x>X0,   t>0
            u(x,0+)   = x*sin(3*x)/6
            u(X0,t)   = (X0+t)*sin(3*(X0+t))/6

    whose analytical solution is

            u(x,t) = (x+t)*sin(3*(x+t))/6

    and Laplace Transform is

                     (s*x+1)*sin(3*x)+3*x*cos(3*x)   s*cos(3*x)-3*sin(3*x)
            U(x,s) = ----------------------------- + ---------------------
                               6*(s^2+9)                   (s^2+9)^2

    where s = +/- 3*i are double poles.
*/

#include <complex.h>


unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0)
/* Abscissa of convergence and singularities of F(s) with their polar multiplicities */
{
    unsigned int NsingsTOT = 2; /* total number of singularities s_j */
    *Nsings = 1;                /* number of singularities s_j with Im(s_j) ge 0 */
    *sigma0 = 0;                /* abscissa of convergence */

    SINGS[0]=0.0+3.0*I;   MULT[0]=2;
    SINGS[1]=0.0-3.0*I;   MULT[1]=2;

    return NsingsTOT;
}
