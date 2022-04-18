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

#include <math.h>

double ILTfun2(double x, double t)
/* Inverse Laplace Transform function u(x,t) */
{   return sin(t)*(1.0-1.0/exp(x));
/*    return sin(t)*(1.0-exp(-x)); */
}
