/**     EXAMPLE 3a

    Application of Talbot's method to solve
    the following PDE problem

        u_t (x,t) = u_xx (x,t),		x>0, t>0
		  u(x,0+) = x*(x-1)
		  u(0,t)  = 2*t
		  u_x(0,t)= -1

	in [0,X1]x[T0,T1].

    The analytical solution is

            u(x,t) = x*(x-1) + 2*t

    and Laplace Transform is

                      2    x*(x-1)
            U(x,s) = --- + -------
                     s^2      s

    where s = 0 is a double pole.
*/


double ILTfun2(double x, double t)
/* Inverse Laplace Transform function u(x,t) */
{
    return  x*(x-1) + 2*t;
}
