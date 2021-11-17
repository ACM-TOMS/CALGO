#include "../../../TalbotSuiteDE/COM/main_include.h"
#include "../../../TalbotSuiteDE/COM/COM_Talbot_pack.h"
#include "../../../TalbotSuiteDE/COM_DE/COM_Talbot_pack_DE.h" /* for the correction to NOPTS                     */
#include "../../../TalbotSuiteDE/FUN_DE/SEQ_Talbot_pack_DE.h"

#include "mex.h"

#include <time.h>               /* for timing                         */

#if defined _WIN64 || defined _WIN32
    #include <windows.h>        /* for Windows high resolution timing */
#endif


/** PROTOTYPES **/

/* Abscissa of convergence and singularities of F(s) with their polar multiplicities: required */
unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0);

double complex *SEQ_LTsamples_ode (unsigned int NXval, double X[], unsigned int NOPTS, double complex S[], double tol);

double *linspace (unsigned int Nv, double Vmin, double Vmax);


/** Porting mex-function [SEQ_PDE_skill_TIMES.c]
    IN MATLAB call it as:
        [PARtime,LTStime,SUMtime,TOTtime] = SEQ_PDE_skill_TIMES (jFUN,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
         1       2       3       4                               1    2     3    4    5     6    7    8     **/
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    unsigned int h, k; /* for loop index */

    /* CHECK FOR THE PROPER NUMBER OF INPUT/OUTPUT ARGUMENTS */
    int NARGin = 8, NARGout = 4; /* number of input and output parameters in MATLAB call */
    if (nrhs != NARGin)
        mexErrMsgTxt("Wrong number of input parameters to SEQ_PDE_skill_TIMES mex function.");

    if (nlhs != NARGout)
        mexErrMsgTxt("Wrong number of output parameters to SEQ_PDE_skill_TIMES mex function.");


    /* EXTRACT INPUT PARAMETERS  */
    unsigned int jFUN  = (unsigned int)mxGetScalar(prhs[0]); /* switch between SEQ_Talbot1 (1) and SEQ_Talbot2 (2) */
    unsigned int NTval = (unsigned int)mxGetScalar(prhs[1]);
    double       Tmin  = mxGetScalar(prhs[2]),
                 Tmax  = mxGetScalar(prhs[3]);
    unsigned int NXval = (unsigned int)mxGetScalar(prhs[4]);
    double       Xmin  = mxGetScalar(prhs[5]),
                 Xmax  = mxGetScalar(prhs[6]);
    double       tol   = mxGetScalar(prhs[7]);

    char *STR = "\n\tLT samples computed by the user defined function: SEQ_LTsamples_ode() [the MATLAB ode45.m/ode113.m function is used]\n";

    /* TALBOT SUITE DE FUNCTION */
    char *Talbot_fun[] = {"SEQ_TalbotSUM1_DE",  /* modified Talbot's method  */
                          "SEQ_TalbotSUM2_DE"}; /* classical Talbot's method */


    /* INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION */
    double sigma0;            /* abscissa of convergence                       */
    unsigned int NsingsTOT;   /* total number of singularities s_j             */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    NsingsTOT = LTsings2(&Nsings, SINGS, MULT, &sigma0);

    /* INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM */
    double *t = linspace(NTval,Tmin,Tmax);

    /* INPUT DATA ABOUT ODE, COMING FROM APPLYING THE LAPLACE METHOD TO PDE */
    double *x = linspace(NXval,Xmin,Xmax);


    /* OUTPUT DATA */
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot, *IFAIL; /* global and local error flags */
    NUMft = (double *)malloc(NXval*NTval*sizeof(double));
    if (NUMft == NULL)
        {mexErrMsgTxt("\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n");}

    IFAIL = (int *)malloc(NXval*NTval*sizeof(int));
    if (IFAIL == NULL)
        {mexErrMsgTxt("\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n");}


	/* VARIABLES OF THIS FUNCTION FOR TIMING */
    double PARtime, LTStime, SUMtime, TOTtime;
    double TOT_TIME;
    unsigned int NTIMES = 2, Nt;
    #ifdef __unix__
        mexPrintf("\nThis is Unix: %u bit C library\n", (unsigned int)sizeof(void*)*8 );
        struct timespec ts1, ts2; double T0, T1;
    #endif

    #if defined(_WIN32) || defined(_WIN64)
        mexPrintf("\nThis is Windows: %u bit C library\n", (unsigned int)sizeof(void*)*8);
        LARGE_INTEGER ticksPerSecond, TICKS1, TICKS2;
        QueryPerformanceFrequency(&ticksPerSecond); /* processor clock frequency */
        mexPrintf("\tCLOCK FREQUENCY by QueryPerformanceFrequency(): %g GHz\n\n",(double)ticksPerSecond.QuadPart/1e6);
    #endif


    /* TALBOT SUITE DE: skill level  */
    double CONLAM, CONSIG, CONNU; unsigned int NOPTS;
    double complex *S, *FS;
    IFAIL_tot = 0; /* initialize */
    int IFAIL_col = 0;
    unsigned int minNOPTS, maxNOPTS;

    switch (jFUN)
    {
    case 1 : /* modified Talbot's method */

        PARtime = 0.0;
        LTStime = 0.0;
        SUMtime = 0.0;
        for (Nt=0; Nt<NTIMES; Nt++) /* repeatitions */
        {

            /* 1.1) Talbot's parameters at the middle point of [Tmin, Tmax] */
            #ifdef __unix__
                clock_gettime(CLOCK_REALTIME, &ts1);
            #endif
            #if defined(_WIN64) || defined(_WIN32)
                QueryPerformanceCounter(&TICKS1);
            #endif
            /* ************************************************************************************** */
            COM_TalbotPAR (sigma0,(Tmin + Tmax)/2,tol,Nsings,SINGS,MULT,&CONLAM,&CONSIG,&CONNU,&NOPTS);
            NOPTS = COM_TalbotNcorr(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4);
            /* ************************************************************************************** */
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
            PARtime += TOT_TIME;


            /* 1.2) Compute sample points on Talbot's contour and the matrix of
                    Laplace Transform samples FS(i,j) on Talbot's contour */
            #ifdef __unix__
                clock_gettime(CLOCK_REALTIME, &ts1);
            #endif
            #if defined(_WIN64) || defined(_WIN32)
                QueryPerformanceCounter(&TICKS1);
            #endif
            /* *********************************************************************************************** */
            S = (double complex*)malloc(NOPTS*sizeof(double complex));
            if ( S == NULL )
                mexErrMsgTxt("\n***   ERROR in SEQ_PDE_skill_TIMES: DYNAMIC ALLOCATION OF S IS FAILED. ***\n");

            double thetak,   piN = 4.0*atan(1.0)/NOPTS; /* pi over NOPTS (step) */
            unsigned int k;
            for (k=1; k<NOPTS; k++)
            {   thetak = k*piN;
                S[k] = CONSIG + CONLAM*thetak/tan(thetak) + I*CONLAM*CONNU*thetak; /* in math.h cot(x) doesn't exist */
            }
            S[0] = CONSIG + CONLAM;
            FS = SEQ_LTsamples_ode(NXval,x,NOPTS,S,tol);
            free(S);
            /* *********************************************************************************************** */
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
            LTStime += TOT_TIME;


            /* 1.3) Compute, for all x and for all t, u(x,t) with the same LT samples */
            #ifdef __unix__
                clock_gettime(CLOCK_REALTIME, &ts1);
            #endif
            #if defined(_WIN64) || defined(_WIN32)
                QueryPerformanceCounter(&TICKS1);
            #endif
            /* ************************************************************************************ */
            IFAIL_tot = SEQ_TalbotSUM1_DE (CONLAM,CONSIG,CONNU,NOPTS,NXval,FS,NTval,t,NUMft,IFAIL);
            free(FS);
            /* ************************************************************************************ */
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
            SUMtime += TOT_TIME;
        }

        PARtime /= NTIMES; /* mean time for parameters */
        LTStime /= NTIMES; /* mean time for LT samples */
        SUMtime /= NTIMES; /* mean time for summations */
        TOTtime = PARtime + LTStime + SUMtime;
        break;


    case 2 : /* classical Talbot's method for DE */

        PARtime = 0.0;
        LTStime = 0.0;
        SUMtime = 0.0;
        minNOPTS = -1; /* max unsigned int value */
        maxNOPTS =  0;

        for (Nt=0; Nt<NTIMES; Nt++) /* repeatitions */
        {
            unsigned int j;
            for (j=0; j<NTval; j++) /* the computational cost of this loop is not considered */
            {
                /* 2.1) TALBOT'S PARAMETERS AT EACH t[j] */
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts1);
                #endif
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS1);
                #endif
                /* ************************************************************************************** */
                COM_TalbotPAR (sigma0, t[j], tol, Nsings, SINGS, MULT, &CONLAM, &CONSIG, &CONNU, &NOPTS);
                /* ************************************************************************************** */
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
                PARtime += TOT_TIME;

                if (NOPTS < minNOPTS) minNOPTS=NOPTS;
                if (NOPTS > maxNOPTS) maxNOPTS=NOPTS;


                /* 2.2) COMPUTE SAMPLE POINTS ON TALBOT'S CONTOUR AND THE MATRIX
                        OF LAPLACE TRANSFORM SAMPLES FS ON TALBOT'S CONTOUR */
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts1);
                #endif
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS1);
                #endif
                /* *********************************************************************************************** */
                S = (double complex*)malloc(NOPTS*sizeof(double complex));
                if ( S == NULL )
                   mexErrMsgTxt("\n***   ERROR in SEQ_PDE_skill_TIMES: DYNAMIC ALLOCATION OF S IS FAILED. ***\n");

                double thetak,   piN = 4.0*atan(1.0)/NOPTS; /* pi over NOPTS (step) */
                unsigned int k;
                for (k=1; k<NOPTS; k++)
                {   thetak = k*piN;
                    S[k] = CONSIG + CONLAM*thetak/tan(thetak) + I*CONLAM*CONNU*thetak; /* in math.h cot(x) doesn't exist */
                }
                S[0] = CONSIG + CONLAM;
                FS = SEQ_LTsamples_ode(NXval,x,NOPTS,S,tol);
                free(S);
                /* *********************************************************************************************** */
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
                LTStime += TOT_TIME;


                /* 2.3) COMPUTE u(x,t), FOR ALL x AND FOR EACH t[j] */
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts1);
                #endif
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS1);
                #endif
                /* ************************************************************************************ */
                IFAIL_col = SEQ_TalbotSUM2_DE (CONLAM,CONSIG,CONNU,NOPTS,NXval,FS,NTval,t[j],j,NUMft,IFAIL);
                IFAIL_tot = IFAIL_tot || IFAIL_col;
                free(FS);
                /* ************************************************************************************ */
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
                SUMtime += TOT_TIME;
            } /* end-for j<NTval */
        } /* end for Nt<NTIMES */

        PARtime /= NTIMES; /* mean time for parameters */
        LTStime /= NTIMES; /* mean time for LT samples */
        SUMtime /= NTIMES; /* mean time for summations */
        TOTtime = PARtime + LTStime + SUMtime;
        break;


    default :
        mexErrMsgTxt("\n*** No function corresponds to this value ***");
    } /* END switch jFUN */


    /* COPY OUTPUTS TO MATLAB VARIABLES */
    plhs[0] = mxCreateDoubleScalar((double)PARtime);
    plhs[1] = mxCreateDoubleScalar((double)LTStime);
    plhs[2] = mxCreateDoubleScalar((double)SUMtime);
    plhs[3] = mxCreateDoubleScalar((double)TOTtime);


    /* PRINT RESULTS */
    mexPrintf("\nInverting the LT fun at NTval = %4d  values of t in [%5.2f, %5.2f]", NTval,Tmin,Tmax);
    mexPrintf("\n                 and at NXval = %4d  values of x in [%5.2f, %5.2f]", NXval,Xmin,Xmax);
    mexPrintf("\n                 input tolerance: tol = %5.2e\n", tol);
    mexPrintf("%s", STR);
    mexPrintf("\n\t***   RESULTS OF SEQUENTIAL TALBOT SUITE DE [function = %s()]   ***\n", Talbot_fun[jFUN-1]);
    mexPrintf("\n\tExecution times:");
    mexPrintf("\n\tPAR time           = %e     for Talbot's parameters", PARtime);
    mexPrintf("\n\tLT sample time     = %e     for LT samples of U(x,s) on Talbot's contour", LTStime);
    mexPrintf("\n\tSUM time           = %e     for the Talbot-Clenshaw summation", SUMtime);
    if (jFUN == 1)
        mexPrintf("\t[NOPTS = %d]", NOPTS);
    else
        mexPrintf("\t[%d <= NOPTS <= %d]", minNOPTS,maxNOPTS);
    mexPrintf("\n\tTotal elapsed time = %e\n", TOTtime);

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

