/**     EXAMPLE 4b

    Application of Talbot's method to solve
    the following PDE problem

            u_tt (x,t) = u_xx (x,t),            0 < x < L,   t>0
            u (x,0+)   = x*sin(3*x)/6
            u_t (x,0+) = sin(3*x)/6 + x*cos(3*x)/2
            u (0,t)    = t*sin(3*t)/6
            u (L,t) =  (L+t)*sin(3*(L+t))/6

    The analytical solution is

            u(x,t) = (x+t)*sin(3*(x+t))/6

    and Laplace Transform is

            U(x,s) = [sin(3*x)+3*x*cos(3*x)+s*x*sin(3*x)]/[6*(s^2+9)]
                    - [3*sin(3*x)-s*cos(3*x)]/(s^2+9)^2

    where s = +/- 3i are double poles.
*/

#include <complex.h>
#include <math.h>


double complex LTfun2(double x, double complex s)
/* Laplace Transform function U(x,s) */
{
    double XXX = 3*x,  SIN3 = sin(XXX), COS3 = cos(XXX);
    double complex S29 = s*s+9;

    return (double complex)(  ( (s*x+1.0)*SIN3 + XXX*COS3 )/(6*S29)
                            + ( s*COS3 - 3*SIN3 )/(S29*S29)  );
}
