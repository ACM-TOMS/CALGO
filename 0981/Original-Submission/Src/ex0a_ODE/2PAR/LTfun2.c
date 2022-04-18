#include <complex.h>

/* LT fun F(s) */
double complex LTfun(double complex s);

/* dummy LT fun F(x,s) for this example */
double complex LTfun2(double x, double complex s)
{   return (*LTfun)( s );
}
