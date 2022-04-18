/** Example 2.2.6 (pag.100)
    Dean G. Duffy - Transform Methods for solving partial differential equations.
    Chapman & Hall/CRC, 2004
 */

#include <math.h>
#include <complex.h>

double complex LTfun(double complex s)
/* Laplace Transform function U(eta,s) */
{
    /** global variable for eta in U(eta,s) **/
    extern double eta;

    double complex LapVal = s*s+9;
    return exp(-eta)*s/(LapVal*LapVal)*cexp(-eta*csqrt(s+1.0));
}
