/********************        SEQ_main_TIMES.c        **********************
 *                                                                        *
 *                                                                        *
 *                        SEQUENTIAL TALBOT SUITE DE                      *
 *                                                                        *
 *                       DRIVER PROGRAM FOR EXAMPLE 5                     *
 *                                TIME TEST                               *
 *                                                                        *
 *                                                                        *
 *                       AUTHOR: Mariarosaria Rizzardi                    *
 *                                                                        *
 *                  mariarosaria.rizzardi@uniparthenope.it                *
 *                                                                        *
 *                  DiST - Dept. of Science and Technology                *
 *                  "Parthenope" University, Naples (Italy)               *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 * REFERENCES                                                             *
 * ==========                                                             *
 * M. RIZZARDI: "Algorithm xxx: TALBOT SUITE DE: APPLICATION OF MODIFIED  *
 *                              TALBOT'S METHOD TO SOLVE DIFFERENTIAL     *
 *                              PROBLEMS".                                *
 *                                 ACM TRANS. MATH. SOFTWARE, VOL. xx,    *
 *                                 NO. x, month year, (PP ##).            *
 *                                                                        *
 **************************************************************************/

#include "../../../TalbotSuiteDE/COM/main_include.h"
#include "../../../TalbotSuiteDE/COM/COM_Talbot_pack.h"       /* for Talbot Suite's shared functions             */
#include "../../../TalbotSuiteDE/COM_DE/COM_Talbot_pack_DE.h" /* for the correction to NOPTS                     */
#include "../../../TalbotSuiteDE/FUN_DE/SEQ_Talbot_pack_DE.h" /* for Talbot Suite DE's sequential implementation */
#include "../../problem_string.h"

#include "mex.h"

#include <time.h>               /* for timing                         */

#if defined _WIN64 || defined _WIN32
    #include <windows.h>        /* for Windows high resolution timing */
#endif


/* Abscissa of convergence and singularities of F(s) with their polar multiplicities: required */
unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0);

/* LT samples on Talbot's contour */
double complex *SEQ_LTsamples_blktrd (unsigned int NXYval, double XY[], unsigned int N, double complex S[], double tol);

double *linspace (unsigned int Nv, double Vmin, double Vmax);


/*  Porting mex-function [SEQ_main_TIMES.c]
                            no output params          4 input params */
void mexFunction (int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
/*  To run from MATLAB:
            >> SEQ_main_TIMES(tol,jFUN,NTval,NXYval)
    Input parameters:
            1) tol    : tolerance on error
            2) jFUN   : number of SEQ function to be tested
                    = 1 SEQ_Talbot1 [ modified Talbot's method]
                    = 2 SEQ_Talbot2 [classical Talbot's method]
            3) NTval  : number of t-values
            4) NXYval : number of values for (x,y)
 */
{
    if (nlhs != 0)
        mexErrMsgTxt( "No output arguments required.");


    /*   INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM              */
    double tol;                             /* tolerance            */
    unsigned int jFUN;        /* SELECTION OF FUNCTION TO BE TESTED */
    unsigned int NTval, NXYval;       /* number of t and x,y values */
    if (nrhs == 4)
    {
        tol   = mxGetScalar(prhs[0]);                /* 1st argument to mexFunction */
        jFUN  = (unsigned int)mxGetScalar(prhs[1]);  /* 2nd argument to mexFunction */
        NTval = (unsigned int)mxGetScalar(prhs[2]);  /* 3rd argument to mexFunction */
        NXYval = (unsigned int)mxGetScalar(prhs[3]); /* 4th argument to mexFunction */
    }
    else
    {
        tol    = 1e-8;                               /* default values              */
        jFUN   = 1;
        NTval  = 20;
        NXYval = 20;
    }


    /* endpoints of intervals */
    double Tmin=100.0,   Tmax=500.0;
    double XYmin=0.0,    XYmax=1.0;


    char *STR = "SEQ_LTsamples_blktrd";       /* sequential LT samples */
    char *Talbot_fun[] = {"SEQ_Talbot1_DE",  /* sequential  modified Talbot's method */
                          "SEQ_Talbot2_DE"}; /* sequential classical Talbot's method */

    unsigned int k; /* for-loop index */

    /* INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION */
    double sigma0;            /* abscissa of convergence                       */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    LTsings2(&Nsings, SINGS, MULT, &sigma0);


    char *LTS = "LT samples by solving in MATLAB block tridiagonal systems";
    mexPrintf("%s", str);
    mexPrintf("\n\t%s\n", LTS);
    mexPrintf("\n\t      t in [%g, %g],    x,y in [%g, %g],    tol=%e", Tmin,Tmax,XYmin,XYmax,tol);
    mexPrintf("\n====================================================================================\n\n");


    #ifdef __unix__
        mexPrintf("This is Unix: %u bit C library\n", (unsigned int)sizeof(void*)*8 );
        struct timespec ts1, ts2; double T0, T1;
    #endif

    #if defined(_WIN32) || defined(_WIN64)
        mexPrintf("This is Windows: %u bit C library\n", (unsigned int)sizeof(void*)*8);
        LARGE_INTEGER ticksPerSecond, TICKS1, TICKS2;
        QueryPerformanceFrequency(&ticksPerSecond); /* processor clock frequency */
        mexPrintf("\tCLOCK FREQUENCY by QueryPerformanceFrequency(): %g GHz\n",(double)ticksPerSecond.QuadPart/1e6);
    #endif


    /* COMPUTE THE t VALUES */
    double *t = linspace(NTval,Tmin,Tmax);


    /* ********************************************************************** *
     * ALLOCATE INPUT/OUTPUT PARAMETERS TO MATLAB internalMeshPts.m FUNCTION  */
    mxArray *plhs1[1], *prhs1[3]; /* pointers to left and right side MATLAB parameters */

    /* BUILD INPUT PARAMETERS TO internalMeshPts.m function */
    prhs1[0] = mxCreateDoubleScalar((double)NXYval);
    prhs1[1] = mxCreateDoubleScalar(XYmin);
    prhs1[2] = mxCreateDoubleScalar(XYmax);

    /*  CALL THE MATLAB FUNCTION internalMeshPts:  XY = internalMeshPts( NXYval, XYmin, XYmax )
        BY MEANS OF:
         mexCallMATLAB(nlhs, plhs, nrhs, prhs,   functionName) */
    if ( mexCallMATLAB(  1,  plhs1,   3,  prhs1, "internalMeshPts") )
        mexErrMsgTxt( "Error in mexCallMATLAB: call to internalMeshPts()");

    mxDestroyArray(prhs1[0]); /* deallocate memory */
    mxDestroyArray(prhs1[1]);
    mxDestroyArray(prhs1[2]);
    /* ********************************************************************** */

    unsigned int Nrows = (NXYval-2)*(NXYval-2); /* number of internal mesh points */

    /*  EXTRACT THE OUTPUT PARAMETER XY OF internalMeshPts():
        XY is a col-wise matrix of size (Nrows, 2)
            XY(i,j) = XY[j*Nrows+i])
        XY = [X Y]
     * ********************************************************************** */
    double  *vr = mxGetPr(plhs1[0]);
    double *XY = (double*)malloc( Nrows*2*sizeof(double) );
    for (k=0; k<2*Nrows; k++)
        XY[k] = *(vr+k);

    mxDestroyArray(plhs1[0]); /* deallocate memory */
    /* ********************************************************************** */

    /*  ALLOCATE OUTPUT DATA
        NUMft, IFAIL: output matrices of size of ((NXYval-2)^2, NTval)) */
    double *NUMft;         /* numerical u(x,y,t) array     */
    int IFAIL_tot, *IFAIL; /* global and local error flags */

    /* ALLOCATE ALGORITHM'S OUTPUT ARRAYS */
    NUMft = (double *)malloc(Nrows*NTval*sizeof(double));
    if (NUMft == NULL)
        mexErrMsgTxt("\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n");

    IFAIL = (int *)malloc(Nrows*NTval*sizeof(int));
    if (IFAIL == NULL)
        mexErrMsgTxt("\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n");


    /* FOR TIMING [both Windows and Linux] */
    double TOT_TIME, MEAN_TIME;
    unsigned int NTIMES = 5, Nt;

    mexPrintf("\n\t***   RESULTS OF SEQ TALBOT SUITE DE function:   %s()   ***", Talbot_fun[jFUN-1]);
    mexPrintf("\n\t      %d t in [%g, %g],    %d x,y in [%g, %g],    input tolerance: tol = %5.2e", NTval,Tmin,Tmax,NXYval,XYmin,XYmax,tol);
    mexPrintf("\n\t      user defined function to compute the LT samples:  %s()\n", STR);
    mexPrintf("\t      %s Talbot's method\n", (jFUN == 1) ? "modified" : "classical");
    mexPrintf("\n\tMEAN ELAPSED TOTAL TIME = ");
    /* SELECTION OF FUNCTION TO BE TESTED */
    switch (jFUN)
    {
        case 1 : /* sequential modified method of Talbot Suite DE */
            MEAN_TIME = 0;
            for (Nt=0; Nt<NTIMES; Nt++)
            {
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts1);
                #endif
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS1);
                #endif
                        IFAIL_tot = SEQ_Talbot1_DE ( SEQ_LTsamples_blktrd, sigma0,Nrows,XY,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax);
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts2);
                    T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec      */
                    T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
                    TOT_TIME = T0 + T1*1.0e-9;
                #endif
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS2);
                    TOT_TIME = (double)(TICKS2.QuadPart - TICKS1.QuadPart)/(double)ticksPerSecond.QuadPart;
                #endif
                MEAN_TIME += TOT_TIME;
            }
            MEAN_TIME /= NTIMES;
            mexPrintf("%e\n\n\n", MEAN_TIME); /* PRINT RESULTS */
            break;

        case 2 : /* classical Talbot's methos */
            MEAN_TIME = 0;
            for (Nt=0; Nt<NTIMES; Nt++)
            {
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts1);
                #endif
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS1);
                #endif
                        IFAIL_tot = SEQ_Talbot2_DE ( SEQ_LTsamples_blktrd, sigma0,Nrows,XY,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT);
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts2);
                    T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec      */
                    T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
                    TOT_TIME = T0 + T1*1.0e-9;
                #endif
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS2);
                    TOT_TIME = (double)(TICKS2.QuadPart - TICKS1.QuadPart)/(double)ticksPerSecond.QuadPart;
                #endif
                MEAN_TIME += TOT_TIME;
            }
            MEAN_TIME /= NTIMES;
            mexPrintf("%e\n\n\n", MEAN_TIME); /* PRINT RESULTS */
            break;

        default :
            mexPrintf( "\n*** No function corresponds to the value %d ***", jFUN);
            exit(1);
    } /* END switch jFUN */
}



double *linspace(unsigned int Nv, double Vmin, double Vmax)
/** Compute a dynamic double array of equispace values **/
{
    double *v = (double *)malloc(Nv*sizeof(double));
    if (v == NULL)
        mexErrMsgTxt("\n***   ERROR IN linspace(): DYNAMIC ALLOCATION OF THE ARRAY IS FAILED. ***\n");

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
