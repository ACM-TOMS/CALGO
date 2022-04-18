/**      EXAMPLE 0

    Laplace Transform function (LTfun) [F24(s) from test functions]:
        F(s) = F(s) = s/(s^2+9)^2

    2 polar singularities (double poles):
        s_j = {+/-3*i}

    abscissa of convergence:
        sigma0 = 0

    Inverse Laplace Transform function (ILTfun):
        f(t) = t*sin(3*t)/6
 */

#include <complex.h>


double complex LTfun(double complex s)
/* Laplace Transform function F(s) */
{   double complex LapVal = s*s+9;
    return s/(LapVal*LapVal);
}
