#include "../../../TalbotSuiteDE/COM/main_include.h"
#include "../../../TalbotSuiteDE/FUN_DE/SEQ_Talbot_pack_DE.h"

#include "mex.h"


/** PROTOTYPES **/
/*  Inverse LT function: not required
    =================================
    This function has been included to compare the numerical approximations
    to the corresponding true results
 */
double ILTfun2 (double x, double t);
unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0);
double *linspace (unsigned int Nv, double Vmin, double Vmax);
double complex *SEQ_LTsamples_ode45 (unsigned int NXval, double X[], unsigned int NOPTS, double complex S[], double tol);



/** Porting mex-function [rel_err.c]
    In MATLAB call it as:
        RELERR = rel_err (jFUN,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
        1                 1    2     3    4    5     6    7    8    number of I/O parameters
    where
        RELERR is a matrix of size (NXval,NTval).
        Its component RELERR(h,k) is the relative error in u(x(h),t(k)).
 **/
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* CHECK FOR THE PROPER NUMBER OF INPUT/OUTPUT ARGUMENTS */
    int NARGin = 8, NARGout = 1; /* number of input and output parameters in MATLAB call */
    if (nrhs != NARGin)
        mexErrMsgTxt("Wrong number of input parameters.");

    if (nlhs != NARGout)
        mexErrMsgTxt("Wrong number of output parameters.");


    /* EXTRACT INPUT PARAMETERS  */
    unsigned int jFUN  = (unsigned int)mxGetScalar(prhs[0]); /* switch between SEQ_Talbot1 (1) and SEQ_Talbot2 (2) */
    unsigned int NTval = (unsigned int)mxGetScalar(prhs[1]);
    double Tmin = mxGetScalar(prhs[2]),
           Tmax = mxGetScalar(prhs[3]);
    unsigned int NXval = (unsigned int)mxGetScalar(prhs[4]);
    double Xmin = mxGetScalar(prhs[5]),
           Xmax = mxGetScalar(prhs[6]);
    double tol  = mxGetScalar(prhs[7]);


    unsigned int h, k; /* for loop index */

    /* INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM */
    double *t = linspace(NTval,Tmin,Tmax);

    /* INPUT DATA ABOUT ODE, COMING FROM APPLYING THE LAPLACE METHOD TO PDE */
    double *x = linspace(NXval,Xmin,Xmax);


    /* OUTPUT DATA */
    double *NUMft;         /* numerical u(x,t) array       */
    int IFAIL_tot, *IFAIL; /* global and local error flags */
    NUMft = (double *)malloc(NXval*NTval*sizeof(double));
    if (NUMft == NULL)
        mexErrMsgTxt("\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n");

    IFAIL = (int *)malloc(NXval*NTval*sizeof(int));
    if (IFAIL == NULL)
        mexErrMsgTxt("\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n");


    /* INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION */
    double sigma0;            /* abscissa of convergence                       */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    LTsings2(&Nsings, SINGS, MULT, &sigma0);


    /* TALBOT SUITE DE user level functions */
    switch (jFUN)
    {
        case 1 : /* MODIFIED TALBOT'S METHOD */
            IFAIL_tot = SEQ_Talbot1_DE ( SEQ_LTsamples_ode45, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax);
            break;

        case 2 : /* CLASSICAL TALBOT'S METHOD */
            IFAIL_tot = SEQ_Talbot2_DE ( SEQ_LTsamples_ode45, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT);
            break;

        default :
            mexErrMsgTxt("\n*** No function corresponds to this value ***");
    } /* END switch jFUN */


    /** COMPUTE THE Relative Error RELERR(h,k) in u(x(h),t(k)) AND COPY IT TO OUTPUT MATLAB MATRIX.
        MATLAB allocates a matrix by cols and C by rows: since RELERR is a C [row-wise] matrix
        we must change its memory allocation in assigning it to the output parameter P=plhs[0].
        For a matrix of size (Nrow, Ncol) we need to set: P(i,j) = RELERR(i,j)
        where
                    RELERR(i,j) = RELERR(i*Ncol+j)  row-wise allocation
                        P(i,j)  = P(i+j*Nrow)       col-wise allocation
    **/
    plhs[0] = mxCreateDoubleMatrix( (mwSize)NXval, (mwSize)NTval, mxREAL);
    double ft,   *P = mxGetPr(plhs[0]);

    for (h=0; h<NXval; h++)     /* row-index */
        for (k=0; k<NTval; k++) /* col-index */
        {
            ft = ILTfun2( x[h], t[k] );
            *(P+h+k*NXval)    = fabs( ft - NUMft[h*NTval + k] )/fabs(ft);
          /* RELERR[h*NTval+k] = fabs( ft - NUMft[h*NTval + k] )/fabs(ft); */
        }

    free(t); free(x); free(NUMft); free(IFAIL);
}



double *linspace(unsigned int Nv, double Vmin, double Vmax)
/** Compute a dynamic double array of equispace values **/
{
    double *v = (double *)malloc(Nv*sizeof(double));
    if (v == NULL)
        {mexErrMsgTxt("\n***   ERROR: DYNAMIC ALLOCATION IS FAILED. ***\n"); exit(1);}

    double step;
    if (Nv>1)
        /* several values */
        step = (Vmax-Vmin)/(Nv-1);
    else
        /* a single value */
        step = 0.0;

    unsigned int k;
    for (k=0; k<Nv; k++)
        *(v+k) = Vmin + k*step;

    return v;
}

