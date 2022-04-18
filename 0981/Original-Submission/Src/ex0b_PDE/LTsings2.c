/** Example 2.2.6 (pag.100)
    Dean G. Duffy - Transform Methods for solving partial differential equations.
    Chapman & Hall/CRC, 2004
 */

#include <complex.h>


unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0)
/* Abscissa of convergence and singularities of F(s) with their polar multiplicities */
{   unsigned int NsingsTOT = 3; /* total number of singularities s_j */
    *Nsings = 2;                /* number of singularities s_j with Im(s_j) ge 0 */
    *sigma0 = 0;                /* abscissa of convergence */

    SINGS[0]=-1.0+0.0*I;    MULT[0]=0; /* branch point */
    SINGS[1]= 0.0+3.0*I;    MULT[1]=2; /* pole */
    SINGS[2]= 0.0-3.0*I;    MULT[2]=2;

    return NsingsTOT;
}
