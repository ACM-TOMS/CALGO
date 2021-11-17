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
#include <math.h>


double complex LTfun2(double x, double complex s)
/* Laplace Transform function U(x,s) */
{
    return (double complex)(1.0 - exp(-x))/(s*s+1);
}
