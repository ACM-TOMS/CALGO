/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR [OpenMP version]
 **/

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#ifdef _OPENMP
    #include <omp.h>
#endif


double complex LTfun2(double x, double complex s);


double complex *OMP_LTsamples_fun (unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS)
/**************************************************************************
    Return the matrix of LT samples computed at (Xval[h], S[k])
    which has to be used row by row in the summation step
**************************************************************************/
{
    /* ALLOCATE OUTPUT ARRAY
       FS must be a matrix of size (NXval,NOPTS)
                FS[h][k] =  U(X[h],S[k])      LaplaceTransform
     */
    double complex *FS = (double complex *)malloc(NXval*NOPTS*sizeof(double complex));
    if (FS == NULL)
    {   fprintf(stderr, "\n***   ERROR IN OMP_LTsamples_fun: DYNAMIC ALLOCATION OF FS IS FAILED. ***\n");
        exit(1);
    }


    /* COMPUTE LT SAMPLES ON TALBOT'S CONTOUR */
    unsigned int j, jS, jX; /* indices */
    unsigned int NTOTloc, STARTloc, ENDloc;
    int          mod, myid;


    unsigned int NTOT = NXval*NOPTS;


    /* PARALLEL SECTION */
    #pragma omp parallel    default   (shared)                  \
                            private   (NTOTloc,STARTloc,ENDloc) \
                            private   (mod,myid,j,jS,jX)        \
                            num_threads (THREADS)
    {
        /* local distribution of indices */
        NTOTloc = NTOT/THREADS;    mod = NTOT%THREADS;
        #ifdef _OPENMP
            myid = omp_get_thread_num();
        #else
            myid = 0;
        #endif
        if (myid < mod)
        {   NTOTloc  = NTOTloc+1;
            STARTloc = myid*NTOTloc;
        }
        else
            STARTloc = myid*NTOTloc + mod;

        ENDloc = STARTloc + NTOTloc - 1;

        /* loop on local values */
        for (j=STARTloc; j<=ENDloc; j++)
        {
            /* COMPUTE LT SAMPLES FS[jX][jS] ON TALBOT'S CONTOUR */
            jX = j/NOPTS; // integer quotient:  jX is the index on x-values
            jS = j%NOPTS; // integer remainder: jS is the index on s-values
            FS[j] = LTfun2(Xval[jX],S[jS]);
        } /* end-loop on local values */

    } /* end parallel section */

    return FS;
}


