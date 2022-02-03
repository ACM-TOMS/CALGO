/********************        OMP_main_ACCURACY.c       ********************
 *                                                                        *
 *                                                                        *
 *                  OpenMP-BASED PARALLEL TALBOT SUITE DE                 *
 *                                                                        *
 *                       DRIVER PROGRAM FOR EXAMPLE 2                     *
 *                              ACCURACY TEST                             *
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
#include "../../../TalbotSuiteDE/FUN_DE/OMP_Talbot_pack_DE.h" /* for Talbot Suite DE's OpenMP implementation */
#include "../../problem_string.h"

#include "mex.h"

#ifdef _OPENMP
    #include <omp.h>
#endif


/*  Inverse LT function u(x,t): not required
    ========================================
    This function has been included to compare the numerical approximations
    to the corresponding true results
 */
double ILTfun2 (double x, double t); /* u(x,t) */

/* Abscissa of convergence and singularities of F(s) with their polar multiplicities: required */
unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0);

/* LT samples on Talbot's contour */
double complex *OMP_LTsamples_ode45 (unsigned int NXval, double Xval[], unsigned int N, double complex S[], double tol, int THREADS);

double *abs_err (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol, int thrds1, int thrds2);
void display_err(char *STR, char FUNtype, unsigned int NXval, unsigned int NTval, double *ERROR, int thrds1, int thrds2);
double *linspace (unsigned int Nv, double Vmin, double Vmax);


/* Porting mex-function [OMP_main_ACCURACY.c]
                            no output params          1 input param */
void mexFunction (int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
/*  THIS PROGRAM REQUIRES THE PARALLEL COMPUTING TOOLBOX (PCT)
    To run from MATLAB:
            >> OMP_main_ACCURACY(tol, thrds1)
    or
            >> OMP_main_ACCURACY(tol, thrds1, threds2)     [for nested parallelism]

    Input parameter:
            1) tol          : (double) tolerance on error
            2) thrds1,thrds2: (int) number of parallel threads (outer and inner respectively)
 */
{
    if (nlhs != 0)
        mexErrMsgTxt( "No output arguments required.");

    double tol;                             /* tolerance            */
    int thrds1, thrds2;                     /* number of OpenMP threads */
    switch (nrhs) {
        case 3 :
            tol = mxGetScalar(prhs[0]);          /* 1st argument to mexFunction */
            thrds1 = (int)mxGetScalar(prhs[1]);  /* 2nd argument to mexFunction */
            thrds2 = (int)mxGetScalar(prhs[2]);  /* 3rd argument to mexFunction */
            break;

        case 2 :
            tol = mxGetScalar(prhs[0]);          /* 1st argument to mexFunction */
            thrds1 = (int)mxGetScalar(prhs[1]);  /* 2nd argument to mexFunction */
            thrds2 = 1;
            break;

        default :
            tol = 1e-8;                         /* default values */
            thrds1 = thrds2 = 1;
    }


    /*  SET THE NUMBER OF OpenMP PARALLEL THREADS */
    #ifdef _OPENMP
        omp_set_dynamic(0); /* static scheduling */
        omp_set_num_threads(thrds1*thrds2);
    #endif

    /* ********************************************************************** *
        SET THE NUMBER OF MATLAB WORKERS
        It is required the MATLAB Parallel Computing Toolbox
        and the number of workers must agree with the local cluster */
    mxArray *plhs1[1], *prhs1[1]; /* pointers to left and right MATLAB params */
    prhs1[0] = mxCreateDoubleScalar((double)thrds1*thrds2); /* 1st input param */
    if ( mexCallMATLAB(  0,  plhs1,   1,  prhs1, "parpool") )
        mexErrMsgTxt( "No parpool available.");
    mxDestroyArray(prhs1[0]); /* deallocate memory */
    /* ********************************************************************** */


    /*   INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM                 */
    unsigned int NXval=9;           /* number of x values              */
    double Xmin= 0.0,  Xmax= 5.0;   /* endpoints of the interval for x */

    unsigned int NTval=5;           /* number of t values              */
    double Tmin=100, Tmax=500;      /* endpoints of the interval for t */


    char *LTS = "LT samples by solving ODE problems by means of MATLAB ode45.m + PCT";
    mexPrintf("%s", str);
    mexPrintf("\n\t   Ex. 2: output from ./2PAR/LTS3_mex/OMP_main_ACCURACY.c\n");
    mexPrintf("\n\t%s\n", LTS);
    mexPrintf("\n\t      %d t in [%g, %g],    %d x in [%g, %g],    tol=%e", NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
    mexPrintf("\n====================================================================================\n\n");


    /* COMPUTE ERRORS */
    double *ABSERR1, *ABSERR2, *ABSERR3;
    ABSERR1 = abs_err (1,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol,thrds1,thrds2);
    ABSERR2 = abs_err (2,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol,thrds1,thrds2);
    ABSERR3 = abs_err (3,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol,thrds1,thrds2);

    /* DISPLAY ERROR MATRICES */
    display_err("ABS",'1',NXval,NTval,ABSERR1,thrds1,thrds2);
    display_err("ABS",'2',NXval,NTval,ABSERR2,thrds1,thrds2);
    display_err("ABS",'3',NXval,NTval,ABSERR3,thrds1,thrds2);

    free(ABSERR1); free(ABSERR2); free(ABSERR3);

    mexPrintf("\n\t      %d t in [%g, %g],    %d x in [%g, %g],    tol = %e\n", NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);


    /* ********************************************************************** *
       DESTROY THE PARALLEL POOL BY CALLING THE MATLAB FUNCTION eval('delete(gcp)') */
    prhs1[0] = mxCreateString("delete(gcp)");
    mexCallMATLAB(  0,  plhs1,  1,  prhs1, "eval");
    mxDestroyArray(prhs1[0]); /* deallocate memory */
    /* ********************************************************************** */
}



double *abs_err (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol, int thrds1, int thrds2)
{
    unsigned int h, k; /* for loop index */

    /* INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM */
    double *t = linspace(NTval,Tmin,Tmax);

    /* INPUT DATA ABOUT ODE, COMING FROM APPLYING THE LAPLACE SUBSTITUTION METHOD TO THE PDE */
    double *x = linspace(NXval,Xmin,Xmax);


    /* OUTPUT DATA */
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot, *IFAIL; /* global and local error flags */
    NUMft = (double *)malloc(NXval*NTval*sizeof(double));
    if (NUMft == NULL)
        {mexErrMsgTxt( "\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n"); exit(1);}

    IFAIL = (int *)malloc(NXval*NTval*sizeof(int));
    if (IFAIL == NULL)
        {mexErrMsgTxt( "\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n"); exit(1);}


    /* INPUT DATA ABOUT LAPLACE TRANSFORM */
    double sigma0;            /* abscissa of convergence                       */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    LTsings2(&Nsings, SINGS, MULT, &sigma0);


    /* TALBOT SUITE DE user level functions  */
    switch (jFUN)
    {
        case 1 : /* modified Talbot's method coarse grain parallelism */
            IFAIL_tot = OMP_Talbot11_DE ( OMP_LTsamples_ode45, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax,thrds1);
            break;

        case 2 : /* modified Talbot's method fine-grain parallelism */
            IFAIL_tot = OMP_Talbot12_DE ( OMP_LTsamples_ode45, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax,thrds1);
            break;

        case 3 : /* modified Talbot's method hybrid fine/coarse-grain parallelism */
            #ifdef _OPENMP
                printf("\nCHECK: NESTED PARALLELISM is %s\n", omp_get_nested() ? "supported" : "not supported");
                if ( !omp_get_nested() )
                {   omp_set_nested(1); /* true */
                    printf("                          Now it is %s\n", omp_get_nested() ? "supported" : "not supported");
                }
            #endif
            IFAIL_tot = OMP_Talbot13_DE ( OMP_LTsamples_ode45, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax,thrds1,thrds2);
            break;

        default :
            mexPrintf( "\n*** No function corresponds to the value %d ***", jFUN);
            exit(1);

    } /* END switch jFUN */

    if ( IFAIL_tot != 0 )
        mexPrintf("Total error indicator: IFAIL_tot = %d", IFAIL_tot);

    /* ABSOLUTE ERRORS */
    double ft;
    double *ABSERR = (double *)malloc(NXval*NTval*sizeof(double));
    if (ABSERR == NULL)
        {mexErrMsgTxt( "\n***   ERROR: DYNAMIC ALLOCATION OF ABSERR IS FAILED. ***\n"); exit(1);}

    for (h=0; h<NXval; h++)
        for (k=0; k<NTval; k++)
        {   ft = ILTfun2( x[h], t[k] );
            ABSERR[h*NTval + k] = fabs( ft - NUMft[h*NTval + k] );
        }

    free(t); free(x); free(NUMft); free(IFAIL);
    return ABSERR;
}


void display_err(char *STR, char FUNtype, unsigned int NXval, unsigned int NTval, double *ERROR, int thrds1, int thrds2)
/** STR = "ABS" or "REL" or "PS_"
    FUNtype = '1', '2' or '3'
    ERROR: row-wise matrix of size (NXval,NTval)
 */
{
    char ch = '%';
    unsigned int h, k; /* for-loop index */

    mexPrintf("\n%sERR%c = [", STR,FUNtype);
    #ifdef _OPENMP
        switch ( FUNtype )
        {
        case '1' :
            mexPrintf(" %c COARSE-GRAIN OMP PARALLELISM with %d threads for modified Talbot's method\n", ch,thrds1 );
            break;
        case '2' :
            mexPrintf(" %c FINE-GRAIN OMP PARALLELISM with %d threads for modified Talbot's method\n", ch,thrds1 );
            break;
        case '3' :
            mexPrintf(" %c NESTED COARSE/FINE-GRAIN OMP PARALLELISM with %d (outer), %d (inner) threads for modified Talbot's method\n", ch,thrds1,thrds2 );
        }
    #else
        mexPrintf("\n");
    #endif
    for (h=0; h<NXval; h++)
    {   mexPrintf("          ");
        for (k=0; k<NTval; k++)    mexPrintf("  %e", ERROR[h*NTval+k]);
        mexPrintf("\n");
    }
    mexPrintf(  "          ];\n");
}


double *linspace(unsigned int Nv, double Vmin, double Vmax)
/** Compute a dynamic double array of equispace values **/
{
    double *v = (double *)malloc(Nv*sizeof(double));
    if (v == NULL)
        {mexErrMsgTxt( "\n***   ERROR IN linspace(): DYNAMIC ALLOCATION OF THE ARRAY IS FAILED. ***\n"); exit(1);}

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
