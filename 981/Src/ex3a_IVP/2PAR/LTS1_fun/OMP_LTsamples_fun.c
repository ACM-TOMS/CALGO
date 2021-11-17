/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR BY SOLVING THE FOLLOWING ODE PROBLEM

            U" = s*U - x*(x-1)
            U(0,s)  = 2/s^2
            U'(0,s) = -1/s

      The analytical solution of the ODE problem is:

            U(x,s) = 2/s^2 + x*(x-1)/s
 **/

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#ifdef _OPENMP
    #include <omp.h>
#endif


/* LT fun U(x,s) */
double complex LTfun2(double x, double complex s);


double complex *OMP_LTsamples_fun (unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS)
/**************************************************************************
    Return the matrix of LT samples computed at (Xval[h], S[k])
    which has to be used row by row in the summation step
**************************************************************************/
{
    /* ALLOCATE OUTPUT ARRAY
       FS must be a matrix of size (NXval,NOPTS)
                FS[h][k] = LaplaceTransform(X[h],S[k])
     */
    double complex *FS = (double complex *)malloc(NXval*NOPTS*sizeof(double complex));
    if (FS == NULL)
    {   fprintf(stderr, "\n***   ERROR IN OMP_LTsamples_fun: DYNAMIC ALLOCATION OF FS IS FAILED. ***\n");
        exit(1);
    }

    unsigned int j, jS, jX; /* indices */
    unsigned int NTOTloc, STARTloc, ENDloc;
    int          mod, myid;


    /* COMPUTE LT SAMPLES FS[h][k] ON TALBOT'S CONTOUR
            FS(h,k) = U(X(h),S(k))
        FS is a row-wise matrix of size (NXval,NOPTS)
       such that
                FS[h][k]  <--->  FS[j],     j=0,...,NTOT-1
       where
            NTOT = NXval*NOPTS
            h=0,...,NXval-1     ==>     h = j/NOPTS
            k=0,...,NOPTS-1     ==>     k = j%NOPTS
    */
    unsigned int NTOT = NXval*NOPTS;

    #pragma omp parallel default   (shared)                  \
                         private   (NTOTloc,STARTloc,ENDloc) \
                         private   (mod,myid,j,jS,jX)        \
                         num_threads (THREADS)
    {
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

        for (j=STARTloc; j<=ENDloc; j++) /* loop on local values */
        {   /*  COMPUTE LT SAMPLES FS[jX][jS] ON TALBOT'S CONTOUR
                jX is the index on x-values
                jS is the index on s-values */
            jX = j/NOPTS;
            jS = j%NOPTS;
            FS[j] = LTfun2 (Xval[jX],S[jS]);
        } /* end-loop on local values */
    } /* end parallel section */

    return FS;
}
