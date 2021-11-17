/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR
 **/

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>


/* LT fun F(s) */
double complex LTfun(double complex s);

/* dummy LT fun U(x,s) for this example */
double complex LTfun2(double x, double complex s)
{   return (*LTfun)( s );
}


double complex *SEQ_LTsamples_fun (unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol)
/**************************************************************************
    Return the matrix of LT samples computed at (Xval[h], S[k])
    which has to be used row by row in the summation step
**************************************************************************/
{
    /* ALLOCATE OUTPUT ARRAY
       FS must be a matrix of size (NXval,NOPTS)
                FS[h][k] = U(X[h],S[k])      LaplaceTransform
     */
    double complex *FS = (double complex *)malloc(NXval*NOPTS*sizeof(double complex));
    if (FS == NULL)
    {   fprintf(stderr, "\n***   ERROR IN SEQ_LTsamples_fun: DYNAMIC ALLOCATION OF FS IS FAILED. ***\n");
        exit(1);
    }

    /* COMPUTE LT SAMPLES ON TALBOT'S CONTOUR */
    unsigned int h, k;    /* for loop indices */


    /* Compute the [row-wise] matrix FS of the LT samples */
    for ( h=0; h<NXval; h++ )
        for ( k=0; k<NOPTS; k++ )
            FS[h*NOPTS+k] = LTfun2(Xval[h],S[k]); /* FS(h,k) */

    return FS;
}


