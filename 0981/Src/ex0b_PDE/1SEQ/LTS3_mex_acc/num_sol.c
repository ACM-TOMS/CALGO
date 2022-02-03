#include "mex.h"

#include "../../../TalbotSuiteDE/COM/main_include.h"
#include "./SEQ_Talbot_pack.h"
#include "../../../TalbotSuiteDE/FUN_DE/SEQ_Talbot_pack_DE.h"


/** global variable for eta in U(eta,s) [ONLY REQUIRED TO USE Talbot Suite] **/
double eta;


/** PROTOTYPES **/
unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0);
double *linspace (unsigned int Nv, double Vmin, double Vmax);

/* F(s): input parameter used by the Talbot Suite functions */
double complex LTfun(double complex s);

/* U(x,s): input parameter used by the Talbot Suite DE functions */
double complex *SEQ_LTsamples_fun (unsigned int NXval, double X[], unsigned int NOPTS, double complex S[], double tol);



/** Porting mex-function [num_sol.c]
    IN MATLAB call it as:
        uxt = num_sol (jFUN,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
        1                1    2     3    4    5     6    7    8    number of I/O parameters
    where
        uxt is a matrix of size (NXval,NTval).
        Its component u(h,k) refers to the numerical approximation of u(x(h),t(k)).
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
    unsigned int jFUN  = (unsigned int)mxGetScalar(prhs[0]);
    /* jFUN = 1  ==>  SEQ_Talbot1
            = 2  ==>  SEQ_Talbot2
            = 3  ==>  SEQ_Talbot1_DE
            = 4  ==>  SEQ_Talbot2_DE
     */
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
    double *NUMft;                      /* numerical u(x,t) array       */
    NUMft = (double *)malloc(NXval*NTval*sizeof(double));
    if (NUMft == NULL)
        mexErrMsgTxt("\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n");

    int IFAIL_tot=0, IFAIL_loc, *IFAIL; /* global and local error flags */
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
        case 1 : /* MODIFIED TALBOT'S METHOD in Talbot Suite */
            for (h=0; h<NXval; h++)
            {
                eta = x[h];
                IFAIL_loc = SEQ_Talbot1 ( LTfun, sigma0,NTval,t,tol,NUMft+h*NTval,IFAIL+h*NTval,Nsings,SINGS,MULT,Tmin,Tmax);
                IFAIL_tot = IFAIL_tot || IFAIL_loc;
            }
            break;

        case 2 : /* CLASSICAL TALBOT'S METHOD in Talbot Suite */
            for (h=0; h<NXval; h++)
            {
                eta = x[h];
                IFAIL_loc = SEQ_Talbot2 ( LTfun, sigma0,NTval,t,tol,NUMft+h*NTval,IFAIL+h*NTval,Nsings,SINGS,MULT);
                IFAIL_tot = IFAIL_tot || IFAIL_loc;
            }
            break;

        case 3 : /* MODIFIED TALBOT'S METHOD in Talbot Suite DE */
            IFAIL_tot = SEQ_Talbot1_DE ( SEQ_LTsamples_fun, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax);
            break;

        case 4 : /* CLASSICAL TALBOT'S METHOD in Talbot Suite DE */
            IFAIL_tot = SEQ_Talbot2_DE ( SEQ_LTsamples_fun, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT);
            break;

        default :
            mexErrMsgTxt("\n*** No function corresponds to this value ***");
    } /* END switch jFUN */


    /** RETURN THE NUMERICAL SOLUTION u(x(h),t(k)) AND COPY IT TO OUTPUT MATLAB MATRIX.
        MATLAB allocates a matrix by cols and C by rows: since NUMft is a C [row-wise] matrix
        we must change its memory allocation in assigning it to the output parameter P=plhs[0].
        For a matrix of size (Nrow, Ncol) we need to set: P(i,j) = NUMft(i,j)
        where
                    NUMft(i,j) = NUMft(i*Ncol+j)  row-wise allocation
                        P(i,j) = P(i+j*Nrow)      col-wise allocation
    **/
    plhs[0] = mxCreateDoubleMatrix( (mwSize)NXval, (mwSize)NTval, mxREAL);
    double ft,   *P = mxGetPr(plhs[0]);

    for (h=0; h<NXval; h++)     /* row-index */
        for (k=0; k<NTval; k++) /* col-index */
            *(P+h+k*NXval) = NUMft[h*NTval + k]; /* function values */

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

