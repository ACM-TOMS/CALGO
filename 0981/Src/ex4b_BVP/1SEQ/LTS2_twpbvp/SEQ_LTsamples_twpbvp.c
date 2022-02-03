/**     USER DEFINED FUNCTION TO COMPUTE THE LAPLACE TRANSFORM SAMPLES
        ON THE TALBOT CONTOUR BY SOLVING THE FOLLOWING ODE PROBLEM

            U"   = s^2*U - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2
            U(0) = s/(s^2+9)^2
                   s*cos(3*L)-3*sin(3*L)   (s*L+1)*sin(3*L)+3*L*cos(3*L)
            U(L) = --------------------- + -----------------------------
                         (s^2+9)^2                    6*(s^2+9)

      The analytical solution is:
                     s*cos(3*x)-3*sin(3*x)   (s*x+1)*sin(3*x)+3*x*cos(3*x)
            U(x,s) = --------------------- + -----------------------------
                           (s^2+9)^2                   6*(s^2+9)

        THE ODE PROBLEM IS SOLVED BY twpbvp.f
**/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>


/*  FORTRAN common blocks (THEIR NAMES FOLLOWED BY AN UNDERSCORE '_')
*/
extern struct /* problem's parameters */
{
    double complex s;
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
                int *nfxpnt, double fixpnt[], int *ntol, int ltol[], double vtol[],
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



double complex *SEQ_LTsamples_twpbvp (unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol)
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
    {   fprintf(stderr, "\n***   ERROR IN SEQ_LTsamples_twpbvp: DYNAMIC ALLOCATION OF FS IS FAILED. ***\n");
        exit(1);
    }


    /** *************  LOCAL WORKING VARIABLES AND PARAMETERS TO twpbvp.f  ************** */
    /* WORKING AREA DIMENSIONS */

    /*  ==>  nmax = 86 */
/*
    static int lwrkfl = 10000,
               lwrkin =  6000,
               MAXxx  =  1000;
*/

    /* Example 3b  ==>  nmax = 129: for NXval=120 */
/*
    static int lwrkfl = 15000,
               lwrkin =  9000,
               MAXxx  =  1000;
*/

    /* Example 4b  ==>  nmax = 259 */
    static int lwrkfl = 30000,
               lwrkin = 18000,
               MAXxx  =  1000;


    if ((int) NXval > MAXxx) /* check maximum size */
    {   fprintf(stderr, "\n***   ERROR IN SEQ_LTsamples_twpbvp: NXval > maximum size = 1000.   ***\n");
        exit(-1);
    }

    /* FORTRAN COMMON BLOCKS INITIALIZATION */
    algprs_.nminit = NXval;   /*  */
    algprs_.pdebug =  0;      /* true or false (for debugging output) */
    algprs_.iprint = -1;      /* Manage output inside twpbvp.f: 0 (default), 1 (more), -1 (none) */

    counts_.nfcall = 0;       /* Number of function evaluations [it may be removed!] */
    probs_.L = Xval[NXval-1]; /* Right end point */

    /* PARAMETERS TO twpbvp.f */
    int     nudim = 4,         /* (I par) The declared row dimension of the array U. nudim >= ncomp */
            ncomp = 4,         /* (I par) Number of eqs and number of components of U at each mesh point */
            nlbc = 2,          /* (I par) Number of boundary conditions at the left endpoint (nlbc <= ncomp) */
            ntol = 4,          /* (I par) Number of tolerances */
            ltol[4],           /* (I par) Array of size ntol. For each i, ltol(i) gives the index of the component
                                          of the computed solution U controlled by the ith tolerance vtol(i) */
            iflbvp,            /* (O par) Output error flag
                                          =  0  success
                                          = -1  one of the input parameters is invalid
                                          =  1  the number of mesh points needed for the next iteration would exceed nmax */
            nmsh = NXval,      /* (I/O par) number of mesh points */
            nmax,              /* (O par) Maximum number of mesh points. nmax <= nucol ??? */
            iwrk[lwrkin],      /* working area */
            nfxpnt = nmsh-2;   /* number of "internal fixed mesh points". If nfxpnt = 0,
                                  only the two endpoints are required to appear in every mesh */

    double  aleft  = Xval[0],  /* left end-point: x in [0,L] */
            aright = Xval[NXval-1], /* right end-point L */
            *fixpnt,           /* (I par) dynamically allocated array of size nfxnt */
            vtol[4],           /* (I par) Array of size ntol. Tolerances */
            xx[MAXxx],
            wrk[lwrkfl];       /* working area */


    int     linear = 1, /* linear = .true. or .false. for a linear or not linear problem */
            giveu  = 0, /* giveu  = .false. */
            givmsh = 1; /* givmsh = .false.
                           givmsh = .true.  ==> user defined mesh points */

    /* U(nudim,MAXxx)
       U(ncomp,MAXxx) : col-wise matrix in FORTRAN. In C it is dynamically allocated
                        and stored by cols:   U(i,j) <---> U[j*ncomp + i] = U[j*nudim + i]
                        Each column contains:   U(1,j) = real(U(xx(j),s))
                                                U(2,j) = imag(u(xx(j),s))
                                                U(3,j) = real(U'(xx(j),s))
                                                U(4,j) = imag(U'(xx(j),s))
     */
    double  *U = (double*)malloc(ncomp*MAXxx*sizeof(double));
    if (U == NULL)
    {   fprintf(stderr,"***   ERROR IN SEQ_LTsamples_twpbvp: dynamic allocation of U failed.   ***");
        exit(-1);
    }

    /*  xx:     array of (user-defined) mesh points.
        fixpnt: array of internal user-defined mesh points.
                It agrees with  Xval[1],...,Xval[NXval-2]
    */
    fixpnt = (double*)malloc(nfxpnt*sizeof(double));
    if (fixpnt == NULL)
    {   fprintf(stderr,"*** Error: dynamic allocation of fixpnt failed. ***");
        exit(-1);
    }
    for (k=1; k<nmsh-1; k++)
    {
        fixpnt[k-1] = Xval[k];
        xx[k] = Xval[k];
    }
    xx[0]=aleft;  xx[nmsh-1]=aright;


    /* ltol[]: array of FORTRAN 77 indices (they are > 0) */
    ltol[0]=1;    vtol[0]=tol;       /* tolerance on the function */
    ltol[1]=2;    vtol[1]=tol;
    ltol[2]=3;    vtol[2]=sqrt(tol); /* tolerance on the derivate */
    ltol[3]=4;    vtol[3]=sqrt(tol);


    /* Compute, for each S[k] on Talbot's contour, an entire column of FS */
    for (k=0; k<(int)NOPTS; k++)
    {
        probs_.s = S[k];
        iflbvp = 0;
        twpbvp_ (&ncomp,&nlbc,&aleft,&aright,&nfxpnt,fixpnt, /* fixpnt: pointer to dynamic array */
                 &ntol,&ltol[0],&vtol[0],
                 &linear,&givmsh,&giveu,&nmsh,
                 &xx[0],&nudim,&U[0],&nmax,&lwrkfl,&wrk[0],&lwrkin,&iwrk[0],
                 fsub_, dfsub_, gsub_, dgsub_,
                 &iflbvp);

        /* Check on error flag */
        if ( iflbvp == 1  &&  nmsh > nmax )
        {
            printf("\tError in twpbvp for k =%3d: nmsh > nmax  ==>  abort!\n", k);
            exit(-1);
        }

        if ( iflbvp == -1 )
        {   fprintf(stderr,"\tError in twpbvp:  k =%3d: iflbvp = -1  ==>  abort!\n", k);
            exit(-1);
        }


        /* Returned bu twpbvp: U(ncomp,nmsh) [col-wise real matrix in FORTRAN and stored by cols in C]
                U(:,j)  <---> U[j*ncomp + 0] + I*U[j*ncomp + 1]
                U'(:,j) <---> U[j*ncomp + 2] + I*U[j*ncomp + 3]
            In order to copy, in a column of FS, the real and imaginary parts of the
            solution U at desired mesh points, we interpolate, if necessary, by means of
            the 3rd degree piece-wise Hermite polynomial
        */
        double step = (xx[nmsh-1] - xx[0])/(nmsh-1);
        double reFS, imFS;
        if (nmsh > (int)NXval) /* if different */
        {
            /* only internal points are interpolated */
            for (j=1; j<(int)NXval-1; j++)
            {
                h = floor( (Xval[j]-xx[0])/step ); /* Xval[j] in [xx[h], xx[h+1]] */
                /*             x0,   x1,     f(x0),         f(x1),             f'(x0),        f'(x1),            x */
                reFS = hermite(xx[h],xx[h+1],U[h*ncomp + 0],U[(h+1)*ncomp + 0],U[h*ncomp + 2],U[(h+1)*ncomp + 2],Xval[j]); /* real */
                imFS = hermite(xx[h],xx[h+1],U[h*ncomp + 1],U[(h+1)*ncomp + 1],U[h*ncomp + 3],U[(h+1)*ncomp + 3],Xval[j]); /* imag */
                FS[j*NOPTS+k] = reFS + I*imFS;
            }

            /* 1st mesh point (j=0) */
            FS[k] = U[0] + I*U[1];

            /* last mesh point (j=NXval-1) */
            FS[(NXval-1)*NOPTS+k] = U[(nmsh-1)*ncomp] + I*U[(nmsh-1)*ncomp+1];

            /* RESET nmsh AND xx[] SINCE THEY HAVE BEEN CHANGED BY twpbvp */
            nmsh = NXval;
            for (j=0; j<(int)NXval; j++)
                xx[j] = Xval[j];
        }
        else /* nmsh == NXval */
            for (h=0; h<(int)NXval; h++)
                FS[h*NOPTS+k] = U[h*ncomp] + I*U[h*ncomp+1];
    }

    free(U);  free(fixpnt);
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

