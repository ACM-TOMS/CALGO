/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR BY SOLVING THE FOLLOWING ODE PROBLEM

            U"   = s*U - x*(x-1),                        0 < x < L
            U(0) = 2/s^2
            U(L) = 2/s^2 + L*(L-1)/s

      The analytical solution is:

            U(x,s) = 2/s^2 + x*(x-1)/s

        THE ODE PROBLEM IS SOLVED BY twpbvp.f
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef _OPENMP
    #include <omp.h>
#endif

/* THIS CONSTANT MUST AGREE WITH THAT ONE IN COMMON BLOCK /probs/ OF testfun.f */
#define maxTHRDS 20


/*  FORTRAN common blocks (THEIR NAMES FOLLOWED BY AN UNDERSCORE '_')
*/
extern struct /* problem's parameters */
{
    double complex st[maxTHRDS];
    double L;
} probs_;

extern struct /* number of function evaluations */
{
    int nfcall;
} counts_;

extern struct /* algorithm's parameters */
{
    int nminit, pdebug, iprint, idum;
    double uval0;
} algprs_;



/*  FORTRAN PROCEDURES (all the params by reference)
    THEIR NAMES FOLLOWED BY AN UNDERSCORE '_'
 */
void twpbvp_   (int *ncomp, int *nlbc, double *aleft, double *aright,
                int *nfxpnt, double fixpnt[], int *ntol, int ltol[], double tol[],
                int *linear, int *givmsh, int *giveu, int *nmsh,
                double xx[], int *nudim, double U[], int *nmax,
                int *lwrkfl, double wrk[], int *lwrkin, int iwrk[],
                void (*fsub_)  (int *ncomp, double *x, double U[], double f[]), /* USER-DEFINED FUNCTIONS */
                void (*dfsub_) (int *ncomp, double *x, double U[], double df[]),
                void (*gsub_)  (int *i, int *ncomp, double U[], double *g),
                void (*dgsub_) (int *i, int *ncomp, double U[], double dg[]),
                int *iflbvp); /* ERROR FLAG INDICATOR */

void fsub_  (int *ncomp, double *x, double U[], double f[]);
void dfsub_ (int *ncomp, double *x, double U[], double df[]);
void gsub_  (int *i, int *ncomp, double U[], double *g);
void dgsub_ (int *i, int *ncomp, double U[], double dg[]);


/* Hermite interpolation function */
double hermite(double x0, double x1, double f0, double f1, double fp0, double fp1, double x);



double complex *OMP_LTsamples_twpbvp (unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS)
/**************************************************************************
    Return the matrix FS of LT samples computed at (Xval[h], S[k])
    by solving ODE problems.
    FS has to be used row by row in the summation step
**************************************************************************/
{
    int h, j, k; /* for loop indices */

    /* ALLOCATE OUTPUT ARRAY
       FS must be a (C row-wise) matrix of size (NXval,NOPTS)
                FS[h][k] = LaplaceTransform(X[h],S[k])
     */
    double complex *FS = (double complex *)malloc(NXval*NOPTS*sizeof(double complex));
    if (FS == NULL)
    {   fprintf(stderr, "\n***   ERROR IN OMP_LTsamples_twpbvp: DYNAMIC ALLOCATION OF FS IS FAILED. ***\n");
        exit(1);
    }


    /** *************  LOCAL WORKING VARIABLES AND PARAMETERS TO twpbvp.f  ************** */
    /* WORKING AREA DIMENSIONS */

    /*  ==>  nmax = 86 */
/*
    int lwrkfl = 10000,
        lwrkin =  6000,
        MAXxx  =  1000;
*/

    /* Example 3b  ==>  nmax = 129: for NXval=120 */
    int lwrkfl = 15000,
        lwrkin =  9000,
        MAXxx  =  1000;

    /* Example 4b  ==>  nmax = 259 */
/*
    int lwrkfl = 30000,
        lwrkin = 18000,
        MAXxx  =  1000;
*/


    if ((int) NXval > MAXxx) /* check maximum size */
    {   fprintf(stderr, "\n***   ERROR IN OMP_LTsamples_twpbvp: NXval > maximum size = %d.   ***\n", MAXxx);
        exit(-1);
    }

    /* FORTRAN COMMON BLOCKS INITIALIZATION */
    algprs_.nminit = NXval;
    algprs_.pdebug =  0;      /* true or false (for debugging output) */
    algprs_.iprint = -1;      /* Manage output inside twpbvp.f: 0 (default), 1 (more), -1 (none) */

    counts_.nfcall = 0;       /* Number of function evaluations [it may be removed!] */
    probs_.L = Xval[NXval-1]; /* Right end point */

    /* PARAMETERS TO twpbvp.f */
    int     nudim = 4,         /* nudim >= ncomp */
            ncomp = 4,         /* number of eqs and number of components in U(x,s) */
            nlbc = 2,          /* number of boundary conditions at the left endpoint (nlbc <= ncomp) */
            ntol = 4,          /* number of tolerances */
            ltol[4],           /* (I par) Array of size ntol. For each i, ltol(i) gives the index of the component
                                          of the computed solution U controlled by the ith tolerance vtol(i) */
            iflbvp,            /* (O par) Output error flag
                                          =  0  success
                                          = -1  one of the input parameters is invalid
                                          =  1  the number of mesh points needed for the next iteration would exceed nmax */
            nmsh = NXval,      /* (I/O par) number of mesh points */
            nmax = MAXxx,      /* (O par) Maximum number of mesh points. nmax <= MAXxx */
            nfxpnt = nmsh-2;   /* number of "internal fixed mesh points". If nfxpnt = 0,
                                  only the two endpoints are required to appear in every mesh */

    double  aleft  = Xval[0],  /* left  end-point: x in [0,L] */
            aright = Xval[NXval-1], /* right end-point L */
            *fixpnt,           /* (I par) dynamically allocated array of size nfxpnt */
            vtol[4];           /* (I par) Array of size ntol. Tolerances */

    int     linear = 1, /* linear = .true. or .false. for a linear or not linear problem */
            giveu  = 0, /* giveu  = .false. */
            givmsh = 1; /* givmsh = .false.
                           givmsh = .true.  ==> user defined mesh points */


    /* U is an (O par) array of pointer, where each component U[k] is a pointer to a matrix.
       U[k] --> U(ncomp,MAXxx) : col-wise matrix in FORTRAN. In C it is dynamically allocated
                        and stored by cols:   U(i,j) <---> U[j*ncomp + i] = U[j*nudim + i]
                        Each column contains:   U(1,j) = real(U(xx(j),s))
                                                U(2,j) = imag(U(xx(j),s))
                                                U(3,j) = real(U'(xx(j),s))
                                                U(4,j) = imag(U'(xx(j),s))
       For thread safe, each thread makes use of its own component.
     */
    double **U = (double**)malloc(THREADS*sizeof(double*));
    if (U == NULL)
    {   fprintf(stderr,"***   ERROR IN OMP_LTsamples_twpbvp: dynamic allocation of U failed.   ***\n");
        exit(-1);
    }

    /* THE SAME HOLDS FOR xx, wrk, iwrk ARRAYS */
    double **xx = (double**)malloc(THREADS*sizeof(double*));
    if (xx == NULL)
    {   fprintf(stderr,"***   ERROR IN OMP_LTsamples_twpbvp: dynamic allocation of xx failed.   ***\n");
        exit(-1);
    }

    double **wrk = (double**)malloc(THREADS*sizeof(double*));
    if (wrk == NULL)
    {   fprintf(stderr,"***   ERROR IN OMP_LTsamples_twpbvp: dynamic allocation of wrk failed.   ***\n");
        exit(-1);
    }

    int **iwrk = (int**)malloc(THREADS*sizeof(int*));
    if (iwrk == NULL)
    {   fprintf(stderr,"***   ERROR IN OMP_LTsamples_twpbvp: dynamic allocation of iwrk failed.   ***\n");
        exit(-1);
    }

    for (k=0; k<THREADS; k++)
    {
        U[k] = (double*)malloc(ncomp*MAXxx*sizeof(double));
        if ( U[k] == NULL )
        {   fprintf(stderr,"***   ERROR IN OMP_LTsamples_twpbvp: dynamic allocation of U[%d] failed.   ***\n", k);
            exit(-1);
        }

        xx[k] = (double*)malloc(MAXxx*sizeof(double));
        if ( xx[k] == NULL )
        {   fprintf(stderr,"***   ERROR IN OMP_LTsamples_twpbvp: dynamic allocation of xx[%d] failed.   ***\n", k);
            exit(-1);
        }

        wrk[k] = (double*)calloc(lwrkfl,sizeof(double));
        if ( wrk[k] == NULL )
        {   fprintf(stderr,"***   ERROR IN OMP_LTsamples_twpbvp: dynamic allocation of wrk[%d] failed.   ***\n", k);
            exit(-1);
        }

        iwrk[k] = (int*)calloc(lwrkin,sizeof(int));
        if ( iwrk[k] == NULL )
        {   fprintf(stderr,"***   ERROR IN OMP_LTsamples_twpbvp: dynamic allocation of iwrk[%d] failed.   ***\n", k);
            exit(-1);
        }
    }

    /*  fixpnt: array of internal user-defined mesh points.
                It agrees with  Xval[1],...,Xval[NXval-2]
    */
    fixpnt = (double*)malloc(nfxpnt*sizeof(double));
    if (fixpnt == NULL)
    {   fprintf(stderr,"*** ERROR IN OMP_LTsamples_twpbvp: dynamic allocation of fixpnt failed. ***");
        exit(-1);
    }
    /* parallelize on nmsh */
    #pragma omp parallel for    default   (shared)   \
                                private   (k)        \
                                num_threads (THREADS)
    for (k=1; k<nmsh-1; k++)
        fixpnt[k-1] = Xval[k];


    /* ltol[]: array of FORTRAN 77 indices (they are > 0) */
    ltol[0]=1;    vtol[0]=tol;       /* tolerance on the function */
    ltol[1]=2;    vtol[1]=tol;
    ltol[2]=3;    vtol[2]=sqrt(tol); /* tolerance on the derivate */
    ltol[3]=4;    vtol[3]=sqrt(tol);


    int myid;
    /*  xx:    array of (user-defined) mesh points. */
    #pragma omp parallel for    default   (shared)   \
                                private   (myid,k)   \
                                num_threads (THREADS)
    for (myid=0; myid<THREADS; myid++)
        for (k=0; k<nmsh; k++)
            *(xx[myid] + k) = Xval[k];


    double step, reFS, imFS;
    /* Compute, for each S[k] on Talbot's contour, an entire column of FS */
    #pragma omp parallel for    default   (shared)                \
                                private   (j,k)                   \
                                private   (myid,nmsh,nmax,iflbvp) \
                                private   (h,step,reFS,imFS)      \
                                num_threads (THREADS)
    for (k=0; k<(int)NOPTS; k++)
    {
        myid = omp_get_thread_num();
        probs_.st[myid] = S[k];
        iflbvp = 0;
        twpbvp_ (&ncomp,&nlbc,&aleft,&aright,&nfxpnt,fixpnt, /* fixpnt: pointer to dynamic array */
                 &ntol,&ltol[0],&vtol[0],
                 &linear,&givmsh,&giveu,&nmsh,
                 xx[myid],&nudim,U[myid],&nmax,&lwrkfl,wrk[myid],&lwrkin,iwrk[myid],
                 fsub_, dfsub_, gsub_, dgsub_,
                 &iflbvp);

        /* Check on error flag */
        if ( iflbvp == 1  &&  nmsh > nmax )
        {   fprintf(stderr,"\tError in twpbvp for thrd#%d:  k =%3d: nmsh > nmax  ==>  abort!\n", myid,k);
            exit(-1);
        }

        if ( iflbvp == -1 )
        {   fprintf(stderr,"\tError in twpbvp for thrd#%d:  k =%3d: iflbvp = -1  ==>  abort!\n", myid,k);
            exit(-1);
        }

        /* Returned by twpbvp: U(ncomp,nmsh) [col-wise real matrix in FORTRAN and stored by cols in C]
                U(:,j)  <---> U[j*ncomp + 0] + I*U[j*ncomp + 1]
                U'(:,j) <---> U[j*ncomp + 2] + I*U[j*ncomp + 3]
            In order to copy, in a column of FS, the real and imaginary parts of the
            solution U at desired mesh points, we interpolate, if necessary, by means of
            the 3rd degree piece-wise Hermite polynomial
        */
        step = (*(xx[myid] + nmsh-1) - *xx[myid])/(nmsh-1);
        if (nmsh > (int)NXval) /* if different */
        {
            /* only internal points are interpolated */
            for (j=1; j<(int)NXval-1; j++)
            {
                h = floor( (Xval[j] - *xx[myid])/step ); /* Xval[j] in [xx[h], xx[h+1]] */
                /*             x0,   x1,     f(x0),         f(x1),             f'(x0),        f'(x1),            x */
                reFS = hermite(*(xx[myid]+h),*(xx[myid]+h+1),*(U[myid] + h*ncomp + 0),*(U[myid] + (h+1)*ncomp + 0),*(U[myid] + h*ncomp + 2),*(U[myid] + (h+1)*ncomp + 2),Xval[j]); /* real */
                imFS = hermite(*(xx[myid]+h),*(xx[myid]+h+1),*(U[myid] + h*ncomp + 1),*(U[myid] + (h+1)*ncomp + 1),*(U[myid] + h*ncomp + 3),*(U[myid] + (h+1)*ncomp + 3),Xval[j]); /* imag */
                FS[j*NOPTS + k] = reFS + I*imFS;
            }

            /* 1st mesh point (j=0) */
            FS[k] = *U[myid] + I* *(U[myid] + 1);

            /* last mesh point (j=NXval-1) */
            FS[(NXval-1)*NOPTS + k] = *(U[myid] + (nmsh-1)*ncomp) + I* *(U[myid] + (nmsh-1)*ncomp + 1);

            /* RESET nmsh AND xx[] SINCE THEY HAVE BEEN CHANGED BY twpbvp */
            nmsh = NXval;
            for (j=0; j<(int)NXval; j++)
                *(xx[myid] + j) = Xval[j];
        }
        else /* nmsh == NXval */
            for (h=0; h<(int)NXval; h++)
                FS[h*NOPTS + k] = *(U[myid] + h*ncomp) + I* *(U[myid] + h*ncomp + 1);

    } /* END PARALLEL for (k=0; k<(int)NOPTS; k++) */


    #pragma omp parallel for    default   (shared)   \
                                private   (myid)     \
                                num_threads (THREADS)
    for (myid=0; myid<THREADS; myid++)
    {
        free(U[myid]);
        free(xx[myid]);
        free(wrk[myid]);
        free(iwrk[myid]);
    }

    free(U);  free(xx);  free(fixpnt);
    free(wrk);  free(iwrk);
    /** ********************************************************************************* */

    return FS;
}


double hermite(double x0, double x1, double f0, double f1, double fp0, double fp1, double x)
/** Compute the 3rd degree Hermite interpolant polynomial
        H3(x) : x in [x0,x1]  and f0=f(x0), f1=f(x1), fp0=f'(x0), fp1=f'(x1)
*/
{
    double R0, R1, R02, R12, H3;
    R0 = (x-x0)/(x1-x0);    R02 = R0*R0;
    R1 = (x1-x)/(x1-x0);    R12 = R1*R1;
	H3 = f0*(1 + 2*R0)*R12 + fp0*(x-x0)*R12 + f1*(1 + 2*R1)*R02 + fp1*(x-x1)*R02;
    return H3;
}

