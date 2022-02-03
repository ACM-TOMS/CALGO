/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR BY SOLVING THE FOLLOWING ODE PROBLEM

            U'    = s*U - x*sin(3*x)/6
            U(X0) = [sin(3*X0)+3*X0*cos(3*X0)+s*X0*sin(3*X0)]/[6*(s^2+9)]
                    - [3*sin(3*X0)-s*cos(3*X0)]/(s^2+9)^2

        The analytical solution is:

            U(x,s) = [sin(3*x)+3*x*cos(3*x)+s*x*sin(3*x)]/[6*(s^2+9)]
                    - [3*sin(3*x)-s*cos(3*x)]/(s^2+9)^2

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
                FS[h][k] = U(X[h],S[k])     LaplaceTransform
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
    int neqn = 2; /* neqn: number of equations in each ODE system */
    int iwork[5];
    double x;
    double work[142]; /*double work[100+21*neqn]; */
    double Ux[2];     /*double Ux[neqn]; */
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
    for (k=0; k<NOPTS; k++)
    {
        /* Compute, for each S[k] on Talbot's contour, an entire column of FS */
        x  = Xval[0];                 /* starting point   */
        FS[k] = init_cond(x,S[k],Ux); /* initial conditions*/
        for ( h=1; h<NXval; h++ ) /* row index */
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

            U' = s*U - x*sin(3*x)/6          [ U' = s*U - u(x,0+) ]

    transformed into a real ODE system of two equations.
**/
{
  Up[0] = creal(S)*U[0] - cimag(S)*U[1] - x*sin(3*x)/6.0;
  Up[1] = cimag(S)*U[0] + creal(S)*U[1];
}


double complex init_cond(double x0, double complex S, double U0[])
/** initial conditions for the ODE problem

            U' = s*U - x*sin(3*x)/6          [ U' = s*U - u(x,0+) ]
            U(X0,s) = [sin(3*X0)+3*X0*cos(3*X0)+s*X0*sin(3*X0)]/[6*(s^2+9)]
                     - [3*sin(3*X0)-s*cos(3*X0)]/(s^2+9)^2
**/
{
    double complex FS;

    /* complex initial value */
    double complex S29 = S*S+9;
    double XXX = 3*x0,  SIN3 = sin(XXX),  COS3 = cos(XXX);
    FS = (double complex)((SIN3 + XXX*COS3+S*x0*SIN3)/(6*S29) - (3*SIN3-S*COS3)/(S29*S29));

    U0[0] = creal(FS); /* real initial conditions for col k (1st row) */
    U0[1] = cimag(FS);

    return FS;
}

