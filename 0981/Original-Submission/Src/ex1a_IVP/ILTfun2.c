/**     EXAMPLE 1a

    Application of Talbot's method to solve
    the following PDE problem

            u_t (x,t) = u_x (x,t),  x>X0,   t>0
            u(x,0+) = x
            u(X0,t) = X0 + t

    whose analytical solution is

            u(x,t) = x + t

    and Laplace Transform is

            U(x,s) = x/s + 1/s^2

    where s = 0 is a double pole.
 */

double ILTfun2(double x, double t)
/* Inverse Laplace Transform function u(x,t) */
{
    return x+t;
}
