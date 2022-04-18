/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR BY SOLVING THE FOLLOWING PDE PROBLEM

            U_xx + U_yy -s*U = -[x*(x-1) + y*(y-1)],    0 < x,y < 1

            U(0,y) = U(1,y) = 4/s^2 + y*(y-1)/s
            U(x,0) = U(x,1) = 4/s^2 + x*(x-1)/s


      THE LT SAMPLES ARE GIVEN BY A FUNCTION.

      OPEN MP-BASED IMPLEMENTATION.
 **/

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#ifdef _OPENMP
    #include <omp.h>
#endif


/* LT fun U(x,s) */
double complex LTfun2(double x, double y, double complex s);


double complex *OMP_LTsamples_fun (unsigned int Nrows, double XY[], unsigned int NOPTS, double complex S[], double tol, int THREADS)
/**************************************************************************
    Return the matrix FS of LT samples computed at (Xval[h], S[k])
    by solving ODE problems.
    FS has to be used row by row in the summation step
  *************************************************************************
    XY is a (col-wise) matrix of Nrows rows and two columns, where Nrows = (NXYval-2)^2,
        XY(h,k) = *(XY + k*Nrows + h)
    XY contains the cartesian coordinates of the internal mesh points.
    They have to be passed to the MATLAB function LT_samples().
**************************************************************************/
{
    double complex *FS = (double complex *)malloc(Nrows*NOPTS*sizeof(double complex));
    if ( FS == NULL )
    {   fprintf(stderr, "\n***   ERROR IN OMP_LTsamples_fun: DYNAMIC ALLOCATION OF FS IS FAILED. ***\n");
        exit(1);
    }

    unsigned int j, jS, jXY; /* indices */
    unsigned int NTOTloc, STARTloc, ENDloc;
    int          mod, myid;

    /* COMPUTE LT SAMPLES FS[h][k] ON TALBOT'S CONTOUR
            FS(h,k) = U(X(h),Y(h),S(k))
        FS is a row-wise matrix of size (Nrows,NOPTS)
       such that
                FS[h][k]  <--->  FS[j],     j=0,...,NTOT-1
       where
            NTOT = Nrows*NOPTS
            h=0,...,Nrows-1     ==>     h = j/NOPTS
            k=0,...,NOPTS-1     ==>     k = j%NOPTS
    */
    unsigned int NTOT = Nrows*NOPTS;

    /* PARALLEL SECTION */
    #pragma omp parallel    default   (shared)                  \
                            private   (NTOTloc,STARTloc,ENDloc) \
                            private   (mod,myid,j,jS,jXY)       \
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

        /* loop on local values */
        for (j=STARTloc; j<=ENDloc; j++)
        {   /*  COMPUTE LT SAMPLES FS[jXY][jS] ON TALBOT'S CONTOUR
                jXY is the index on x-values
                jS is the index on s-values */
            jXY = j/NOPTS;
            jS = j%NOPTS;
            /*  X(jXY) = XY(jXY,0)= *(XY + jXY)
                Y(jXY) = XY(jXY,1)= *(XY + Nrows + jXY) */
            FS[j] = LTfun2 (*(XY+jXY), *(XY+Nrows+jXY), S[jS]);
        } /* end-loop on local values */

    } /* end parallel section */

    return FS;
}

