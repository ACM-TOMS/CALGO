/**     EXAMPLE 4a

    Application of Talbot's method to solve
    the following PDE problem

            u_tt (x,t) = u_xx (x,t),            0 < x < L,   t>0
            u (x,0+)   = x*sin(3*x)/6
            u_t (x,0+) = sin(3*x)/6 + x*cos(3*x)/2
            u (0,t)    = t*sin(3*t)/6
            u_x (0,t)  =  sin(3*t)/6 + t*cos(3*t)/2

    The analytical solution is

            u(x,t) = (x+t)*sin(3*(x+t))/6

    and Laplace Transform is

            U(x,s) = [sin(3*x)+3*x*cos(3*x)+s*x*sin(3*x)]/[6*(s^2+9)]
                    - [3*sin(3*x)-s*cos(3*x)]/(s^2+9)^2

    where s = +/- 3i are double poles.
 **/

#include <math.h>


double ILTfun2(double x, double t)
/* Inverse Laplace Transform function u(x,t) */
{
    return (x+t)*sin(3*(x+t))/6.0;
}

