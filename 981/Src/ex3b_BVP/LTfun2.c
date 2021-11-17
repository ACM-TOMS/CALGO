/**     EXAMPLE 3b

    Application of Talbot's method to solve
    the following PDE problem

        u_t (x,t) = u_xx (x,t),		0 < x < L, t > 0
		  u(x,0+) = x*(x-1)
		  u(0,t)  = 2*t
		  u(L,t)  = 2*t + L*(L-1)

    in [0,L]x[T0,T1].

    The analytical solution is

            u(x,t) = x*(x-1) + 2*t

    and Laplace Transform is

                     x*(x-1)    2
            U(x,s) = ------- + ---
                        s      s^2

    where s = 0 is a double pole.
*/

#include <complex.h>
#include <math.h>


double complex LTfun2(double x, double complex s)
/* Laplace Transform function U(x,s) */
{
    return 2.0/(s*s) + x*(x-1)/s;
}
