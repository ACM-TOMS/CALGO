/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR BY SOLVING THE FOLLOWING ODE PROBLEM

            U"    = s^2*U - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2
            U(0)  = s/(s^2+9)^2
            U'(0) = s^2/(s^2+9)^2

      The analytical solution of the ODE problem is:

                     sin(3*x)+3*x*cos(3*x)+s*x*sin(3*x)   3*sin(3*x)-s*cos(3*x)
            U(x,s) = ---------------------------------- - ---------------------
                                  6*(s^2+9)                     (s^2+9)^2

        THE ODE PROBLEM IS SOLVED BY ode.c
 **/

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h> /* for sin(x) */

#include "ode.h" /* for ode.c */

#ifdef _OPENMP
    #include <omp.h>
#endif


/* FUNCTIONS FOR ODE PROBLEM */
void              ODEfun (double x, double U[], double Up[], double complex S);
double complex init_cond (double x0, double complex S, double U0[]);


double complex *OMP_LTsamples_ode (unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS)
/**************************************************************************
    Return the matrix FS of LT samples computed at (Xval[h], S[k])
    by solving ODE problems.
    FS has to be used row by row in the summation step
**************************************************************************/
{
    /* ALLOCATE OUTPUT ARRAY
       FS must be a matrix of size (NXval,NOPTS)
                FS[h][k] = LaplaceTransform(X[h],S[k])
     */
    double complex *FS = (double complex *)malloc(NXval*NOPTS*sizeof(double complex));
    if (FS == NULL)
    {   fprintf(stderr, "\n***   ERROR IN OMP_LTsamples_ode: DYNAMIC ALLOCATION OF FS IS FAILED. ***\n");
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
    int neqn = 4; /* number of equations in each ODE system */
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
        for ( h=1; h<NXval; h++ )
        {
            iflag = 1; /* Initialize the flag indicator */
            ode (ODEfun,S[k],neqn,Ux,&x,Xval[h],tol,tol,&iflag,work,iwork);
            if ( iflag != 2 )
            {   fprintf(stderr, "\n***   ERROR IN OMP_LTsamples_ode: from ode() function iflag = %d. ***\n", iflag);
                exit(1);
            }
            FS[h*NOPTS+k] = Ux[0] + I*Ux[1]; /* FS[h][k] */
        }
    }

    return FS;
}


void ODEfun ( double x, double U[], double Up[], double complex S )
/** complex ODE equation:

            U" = s^2*U - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2
            U(0,s)  = s/(s^2+9)^2
            U'(0,s) = s^2/(s^2+9)^2

    transformed into a real ODE system of four equations:
            Y1+i*Y2 = U1 = U,   Y3+i*Y4 = U2 = U'
    ==>
        Y1' = Y3
        Y2' = Y4
        Y3' = Re[s^2]*Y1 - Im[s^2]*Y2 - (Re[s]*x+1)*sin(3*x)/6 - x*cos(3*x)/2
        Y4' = Im[s^2]*Y1 + Re[s^2]*Y2 - Im[s]*x*sin(3*x)/6
    Real initial conditions
        Y1(0,s) = Re[s/(s^2+9)^2]
        Y2(0,s) = Im[s/(s^2+9)^2]
        Y3(0,s) = Re[s^2/(s^2+9)^2] =
                = Re[s]*Y1(0,s) - Im[s]*Y2(0,s)
        Y4(0,s) = Im[s^2/(s^2+9)^2] =
                = Im[s]*Y1(0,s) + Re[s]*Y2(0,s)
**/
{
    double complex SS = S*S;
    Up[0] = U[2];
    Up[1] = U[3];
    Up[2] = creal(SS)*U[0] - cimag(SS)*U[1] - (creal(S)*x+1)*sin(3*x)/6 - x*cos(3*x)/2;
    Up[3] = cimag(SS)*U[0] + creal(SS)*U[1] -  cimag(S)*x*sin(3*x)/6;
}


double complex init_cond(double x0, double complex S, double U0[])
/** initial conditions for the ODE problem

            U" = s^2*U - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2
            U(0,s)  = s/(s^2+9)^2
            U'(0,s) = s^2/(s^2+9)^2

            Y1+i*Y2 = U1 = U,   Y3+i*Y4 = U2 = U'
    ==>
    Real initial conditions
        Y1(0,s) = Re[s/(s^2+9)^2]
        Y2(0,s) = Im[s/(s^2+9)^2]
        Y3(0,s) = Re[s * s/(s^2+9)^2] =
        Y4(0,s) = Im[s * s/(s^2+9)^2] =
**/
{
    double complex FS;
    /* real initial conditions */
    FS = S/((S*S+9)*(S*S+9));
        U0[0] = creal(FS);
        U0[1] = cimag(FS);
    FS = S*FS;
        U0[2] = creal(FS);
        U0[3] = cimag(FS);

    FS = U0[0] + I*U0[1];
return FS;
}

