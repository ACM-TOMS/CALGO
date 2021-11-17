/********************        SEQ_main_ACCURACY.c       ********************
 *                                                                        *
 *                                                                        *
 *                        SEQUENTIAL TALBOT SUITE DE                      *
 *                                                                        *
 *                       DRIVER PROGRAM FOR EXAMPLE 5                     *
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
#include "../../../TalbotSuiteDE/FUN_DE/SEQ_Talbot_pack_DE.h" /* for Talbot Suite DE's sequential implementation */
#include "../../problem_string.h"

#include "mex.h"


/*  Inverse LT function u(x,t): not required
    ========================================
    This function has been included to compare the numerical approximations
    to the corresponding true results
 */
double ILTfun2 (double x, double y, double t); /* u(x,y,t) */

/* Abscissa of convergence and singularities of F(s) with their polar multiplicities: required */
unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0);

/* LT samples on Talbot's contour */
double complex *SEQ_LTsamples_blktrd (unsigned int NXYval, double XY[], unsigned int N, double complex S[], double tol);

double *rel_err (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXYval, double XYmin, double XYmax, double tol);
double *linspace (unsigned int Nv, double Vmin, double Vmax);


/* Porting mex-function [OMP_main_ACCURACY.c]
                            no output params          input params */
void mexFunction (int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
/*  To call from MATLAB:
            >> SEQ_main_ACCURACY(tol)
    Input parameter:
            1) tol          : (double) tolerance on error
 */
{
    if (nlhs != 0)
        mexErrMsgTxt( "No output arguments required.");


    double tol;                             /* tolerance            */
    if (nrhs == 1)
        tol = mxGetScalar(prhs[0]);          /* 1st argument to mexFunction */
    else
        tol = 1e-8;                         /* default value */


    /*   INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM              */
    unsigned int NXYval=9;       /* number of x,y values            */
    unsigned int  NTval=5;       /* number of t values              */

    /* endpoints */
    double Tmin=100.0,    Tmax=500.0;
    double XYmin=0.0,     XYmax=1.0;


    char *LTS = "LT samples by solving in MATLAB block tridiagonal systems";
    mexPrintf("%s", str);
    mexPrintf("\n\t   Ex. 5: output from ./1SEQ/LTS3_mex_Cmain/SEQ_main_ACCURACY.c");
    mexPrintf("\n\t%s\n", LTS);
    mexPrintf("\n\t      %d t in [%g, %g],    %d x,y in [%g, %g],    tol=%e", NTval,Tmin,Tmax,NXYval,XYmin,XYmax,tol);
    mexPrintf("\n====================================================================================\n\n");


    /* COMPUTE RELATIVE ERRORS */
    double *RELERR1, *RELERR2;
    RELERR1 = rel_err (1,NTval,Tmin,Tmax,NXYval,XYmin,XYmax,tol); /* SEQ_Talbot1_DE */
    RELERR2 = rel_err (2,NTval,Tmin,Tmax,NXYval,XYmin,XYmax,tol); /* SEQ_Talbot2_DE */


    /* DISPLAY ERROR MATRICES */
    unsigned int Nrows = (NXYval-2)*(NXYval-2);
    unsigned int h, k; /* for-loop index */
    char ch = '%';
    mexPrintf("\nRELERR1 = [ %c  Tval(1)    Tval(2)    ...    Tval(%d)\n", ch,NTval);
    for (h=0; h<Nrows; h++)
    {   mexPrintf("%2d        ", h+1);
        for (k=0; k<NTval; k++)    mexPrintf("  %e", RELERR1[h*NTval+k]);
        mexPrintf("    %c XYval(%2d)\n", ch,h+1);
    }
    mexPrintf(  "          ];\n");

    mexPrintf("\nRELERR2 = [ %c  Tval(1)    Tval(2)    ...    Tval(%d)\n", ch,NTval);
    for (h=0; h<Nrows; h++)
    {   mexPrintf("%2d        ", h+1);
        for (k=0; k<NTval; k++)    mexPrintf("  %e", RELERR2[h*NTval+k]);
        mexPrintf("    %c XYval(%2d)\n", ch,h+1);
    }
    mexPrintf(  "          ];\n");
}



double *rel_err (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXYval, double XYmin, double XYmax, double tol)
{
    unsigned int h, k; /* for loop index */

    /* COMPUTE THE t VALUES */
    double *t = linspace(NTval,Tmin,Tmax);


    /* ********************************************************************** *
     * ALLOCATE INPUT/OUTPUT PARAMETERS TO MATLAB internalMeshPts.m FUNCTION  */
    mxArray *plhs[1], *prhs[3]; /* pointers to left and right side MATLAB parameters */

    /* BUILD INPUT PARAMETERS TO internalMeshPts.m function */
    prhs[0] = mxCreateDoubleScalar((double)NXYval);
    prhs[1] = mxCreateDoubleScalar(XYmin);
    prhs[2] = mxCreateDoubleScalar(XYmax);

    /*  CALL THE MATLAB FUNCTION internalMeshPts:  XY = internalMeshPts( NXYval, XYmin, XYmax )
        BY MEANS OF:
        mexCallMATLAB(nlhs, plhs, nrhs, prhs,   functionName) */
        if ( mexCallMATLAB(  1,  plhs,   3,  prhs, "internalMeshPts") )
            mexErrMsgTxt( "Error in mexCallMATLAB: call to internalMeshPts()");

    mxDestroyArray(prhs[0]); /* deallocate memory */
    mxDestroyArray(prhs[1]);
    mxDestroyArray(prhs[2]);
    /* ********************************************************************** */

    unsigned int Nrows = (NXYval-2)*(NXYval-2); /* number of internal mesh points */

    /*  EXTRACT THE OUTPUT PARAMETER XY OF internalMeshPts():
        XY is a col-wise matrix of size (Nrows, 2)
            XY(i,j) = XY[j*Nrows+i])
        XY = [X Y]
    */
    double  *vr = mxGetPr(plhs[0]);
    double *XY = (double*)malloc( Nrows*2*sizeof(double) );
    for (h=0; h<2*Nrows; h++)
        XY[h] = *(vr+h);


    /*  ALLOCATE OUTPUT DATA
        NUMft, IFAIL: output matrices of size of ((NXYval-2)^2, NTval)) */
    double *NUMft;         /* numerical u(x,y,t) array     */
    int IFAIL_tot, *IFAIL; /* global and local error flags */
    NUMft = (double *)malloc(Nrows*NTval*sizeof(double));
    if (NUMft == NULL)
        mexErrMsgTxt( "\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n");

    IFAIL = (int *)malloc(Nrows*NTval*sizeof(int));
    if (IFAIL == NULL)
        mexErrMsgTxt( "\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n");


    /* INPUT DATA ABOUT LAPLACE TRANSFORM */
    double sigma0;            /* abscissa of convergence                       */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    LTsings2(&Nsings, SINGS, MULT, &sigma0);


    /* TALBOT SUITE DE user level functions */
    switch (jFUN)
    {
        case 1 : /* (modified Talbot's method */
            IFAIL_tot = SEQ_Talbot1_DE ( SEQ_LTsamples_blktrd, sigma0,Nrows,XY,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax);
            break;

        case 2 : /* classical Talbot's method */
            IFAIL_tot = SEQ_Talbot2_DE ( SEQ_LTsamples_blktrd, sigma0,Nrows,XY,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT);
            break;

        default :
            mexPrintf( "\n*** No function corresponds to the value %d ***", jFUN);
            exit(1);

    } /* END switch jFUN */

    if ( IFAIL_tot != 0 )
        mexPrintf("Total error indicator: IFAIL_tot = %d", IFAIL_tot);


    /* RELATIVE ERRORS */
    double ft;
    double *RELERR = (double *)malloc(Nrows*NTval*sizeof(double));
    if (RELERR == NULL)
        {mexErrMsgTxt( "\n***   ERROR: DYNAMIC ALLOCATION OF RELERR IS FAILED. ***\n"); exit(1);}

    for (h=0; h<Nrows; h++)
        for (k=0; k<NTval; k++)
        {   ft = ILTfun2( XY[h], XY[Nrows+h], t[k] ); /* col-wise XY(h,j) = XY[j*Nrows+h]) */
            RELERR[h*NTval + k] = fabs( ft - NUMft[h*NTval + k] )/fabs(ft);
        }

    free(t); free(XY); free(NUMft); free(IFAIL);
    return RELERR;
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
