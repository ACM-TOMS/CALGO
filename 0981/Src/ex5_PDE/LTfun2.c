/**     EXAMPLE 5

    Application of Talbot's method to solve
    the following PDE problem

            u_t  = u_xx  + u_yy,            0 < x,y < 1,   t>0
            u(x,y,0+) = x*(x-1) + y*(y-1)
            u(0,y,t) = u(1,y,t) = 4*t + y*(y-1)
            u(x,0,t) = u(x,1,t) = 4*t + x*(x-1)

    The analytical solution is

            u(x,y,t) = 4*t + x*(x-1) + y*(y-1)

    and Laplace Transform is

                      4    x*(x-1) + y*(y-1)
            U(x,s) = --- + -----------------
                     s^2           s

    where s = 0 is a double pole.
**/

#include <complex.h>
#include <math.h>


double complex LTfun2(double x, double y, double complex s)
/* Laplace Transform function U(x,y,s) */
{
    return (double complex)(  4.0 / (s*s) + ( x*(x-1) + y*(y-1) )/s  );
}
