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

#include <math.h>


double ILTfun2(double x, double t)
/* Inverse Laplace Transform function u(x,t) */
{
    return (x+t)*sin(3*(x+t))/6.0;
}

