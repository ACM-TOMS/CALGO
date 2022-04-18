/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR BY SOLVING THE FOLLOWING ODE PROBLEM

                U' = exp(-x)/(s^2+1),        x>0
                U(0,s) = 0

        THE ODE PROBLEM IS SOLVED BY ode.c
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "ode.h" /* for ode.c */

#ifdef _OPENMP
    #include <omp.h>
#endif


/* FUNCTIONS FOR ODE PROBLEM */
void              ODEfun (double x, double U[], double Up[], double complex S);
double complex init_cond (double x0, double complex S, double U0[]);


double complex *OMP_LTsamples_ode (unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS)
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
    {   fprintf(stderr, "\n***   ERROR IN SEQ_LTsamples_ode: DYNAMIC ALLOCATION OF FS IS FAILED. ***\n");
        exit(1);
    }

    unsigned int h, k; /* for loop indices */


    /* *************  LOCAL WORKING VARIABLES FOR ode.c   *************
       Unlike the sequential version, that allocates the working variables
       dynamically in the heap memory:
            double *work = (double*) malloc ( (100+21*neqn)*sizeof(double) );
            double *Ux   = (double*) malloc( neqn*sizeof(double) );
       here, all the working variables are "local" array for thread safety.
    */
    int neqn = 4; /* neqn: number of equations in each ODE system */
    int iwork[5];
    double work[184]; /* double work[100+21*neqn]; */
    double Ux[4];     /* double Ux[neqn] */
    double x;
    /* **************************************************************** */

    int iflag;         /* error indicator from ode.c */


    /* Compute the LT samples in FS, a row-wise matrix of size (NXval,NOPTS),
       such that
                FS[h][k]  <--->  FS[j],     j=0,...,Ntot-1
       where
            Ntot = NXval*NOPTS
            h=0,...,NXval-1     ==>     h = j/NOPTS
            k=0,...,NOPTS-1     ==>     k = j%NOPTS
     */
    #pragma omp parallel for    default   (shared)                    \
                                private   (k,h,Ux,x,iflag,work,iwork) \
                                num_threads (THREADS)
    for (k=0; k<NOPTS; k++) /* parallelized for on column index */
    {
        /* Compute, for each S[k] on Talbot's contour, an entire column of FS */
        x  = Xval[0];               /* starting point   */
        FS[k] = (*init_cond)(x,S[k],Ux); /* intial conditions*/
        for ( h=1; h<NXval; h++ ) /* row index */
        {
            iflag = 1; /* Initialize the flag indicator */
            ode (ODEfun,S[k],neqn,Ux,&x,Xval[h],tol,tol,&iflag,work,iwork);
            if ( iflag != 2 )
                {fprintf(stderr, "\n***   ERROR IN OMP_LTsamples_ode: from ode() function iflag = %d. ***\n", iflag);  exit(1);}
            FS[h*NOPTS+k] = Ux[0] + I*Ux[1]; /* FS[h][k] */
        }
    }

    return FS;
}


void ODEfun ( double x, double U[], double Up[], double complex S )
/** complex ODE equation:

            U" = s*U - x*(x-1)
            U(0,s)  = 2/s^2
            U'(0,s) = -1/s

    transformed into a real ODE system of four equations:
            Y1+i*Y2 = U1 = U,   Y3+i*Y4 = U2 = U'
    ==>
            Y1'+i*Y2' = U1' = U2 = Y3+i*Y4
            Y3'+i*Y4' = U2' = s*U1 - x*(x-1) = s*[Y1+i*Y2] - x*(x-1)

            Y1' = Y3
            Y2' = Y4
            Y3' = Re[s]*Y1 - Im[s]*Y2 - x*(x-1)
            Y4' = Im[s]*Y1 + Re[s]*Y2
**/
{
    Up[0] = U[2];
    Up[1] = U[3];
    Up[2] = creal(S)*U[0] - cimag(S)*U[1] - x*(x-1.0);
    Up[3] = cimag(S)*U[0] + creal(S)*U[1];
}


double complex init_cond(double x0, double complex S, double U0[])
/** initial conditions for the ODE problem

            U" = s*U - x*(x-1)
            U(0,s)  = 2/s^2
            U'(0,s) = -1/s

            Y1+i*Y2 = U1 = U,   Y3+i*Y4 = U2 = U'
    ==>
            Y1(0,s) = Re[2/s^2]
            Y2(0,s) = Im[2/s^2]
            Y3(0,s) = Re[-1/s]
            Y4(0,s) = Im[-1/s]
**/
{
    double complex FS;
    /* real initial conditions */
    FS = 2.0/(S*S);
        U0[0] = creal(FS);
        U0[1] = cimag(FS);
    FS = -1/S;
        U0[2] = creal(FS);
        U0[3] = cimag(FS);
    FS = U0[0] + I*U0[1];
return FS;
}
