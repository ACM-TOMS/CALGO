/********************        SEQ_main_TIMES.c        **********************
 *                                                                        *
 *                                                                        *
 *                        SEQUENTIAL TALBOT SUITE DE                      *
 *                                                                        *
 *                       DRIVER PROGRAM FOR EXAMPLE 2                     *
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
double complex *SEQ_LTsamples_ode45 (unsigned int NXval, double Xval[], unsigned int N, double complex S[], double tol);

double *linspace (unsigned int Nv, double Vmin, double Vmax);
void SEQ_PDE_skill_MAIN (double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol), /* user-defined function for LT samples */
                         unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax,
                         double tol, double *PARtime, double *LTStime, double *SUMtime, double *TOTtime);



/* Porting mex-function [SEQ_main_TIMES.c]
                            no output params          1 input param */
void mexFunction (int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
/*  To call from MATLAB:
            >> SEQ_main_TIMES(1e-10)

    Input parameter:
            tol: tolerance on error
 */
{
    if (nlhs != 0)
        mexErrMsgTxt( "No output arguments required.");

    double tol;                             /* tolerance            */
    if (nrhs == 1)
        tol = mxGetScalar(prhs[0]);         /* 1st argument to mexFunction */
    else
        tol = 1e-8;                         /* default value        */


    /*   INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM              */
    const unsigned int NX=3, NT=3;
    unsigned int NTval[]={5,20,120},   NXval[]={5,20,120};       /* number of x, t values */
    double Tmin=100.0,   Tmax=500.0,   Xmin= 0.0,   Xmax= 5.0; /* endpoints */


    char *LTS = "LT samples by solving ODE problems by means of MATLAB ode45.m";
    mexPrintf("%s", str);
    mexPrintf("\n\t   Ex. 2: output from ./1SEQ/LTS3_mex/SEQ_main_TIMES.c\n");
    mexPrintf("\n\t%s\n", LTS);
    mexPrintf("\n\tt in [%g, %g],    x in [%g, %g],    tol=%e", Tmin,Tmax,Xmin,Xmax,tol);
    mexPrintf("\n====================================================================================\n\n");



    /*  XXXtimeY(h,k): mean execution time T(h,k) for NTval(h), NXval(k) for each step in algorithm.
        XXX: PAR (parameters),
             LTS (LT samples),
             SUM (Talbot's summation)
        Y:   1 (SEQ_Talbot1),
             2 (SEQ_Talbot2)
     */
    double PARtime1[NT][NX], LTStime1[NT][NX], SUMtime1[NT][NX], TOTtime1[NT][NX],
           PARtime2[NT][NX], LTStime2[NT][NX], SUMtime2[NT][NX], TOTtime2[NT][NX];


    unsigned int h, k; /* for-loop index */
    char ch = '%';


    /* COMPUTE THE TIME T(h,k) OF EACH STEP */
    for (h=0; h<NT; h++)        /* Tval[h] */
        for (k=0; k<NX; k++)    /* Xval[k] */
        {
            SEQ_PDE_skill_MAIN (SEQ_LTsamples_ode45,1,NTval[h],Tmin,Tmax,NXval[k],Xmin,Xmax,tol,&PARtime1[h][k],&LTStime1[h][k],&SUMtime1[h][k],&TOTtime1[h][k]);
            SEQ_PDE_skill_MAIN (SEQ_LTsamples_ode45,1,NTval[h],Tmin,Tmax,NXval[k],Xmin,Xmax,tol,&PARtime1[h][k],&LTStime1[h][k],&SUMtime1[h][k],&TOTtime1[h][k]);
            /* chiama 2 volte Talbot1 per ammortizzare l'inizializzazione dei cicli for */
            SEQ_PDE_skill_MAIN (SEQ_LTsamples_ode45,2,NTval[h],Tmin,Tmax,NXval[k],Xmin,Xmax,tol,&PARtime2[h][k],&LTStime2[h][k],&SUMtime2[h][k],&TOTtime2[h][k]);
        }


    /* DISPLAY TIME'S MATRICES */
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nPARtime1 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", PARtime1[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nLTStime1 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", LTStime1[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nSUMtime1 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", SUMtime1[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nTOTtime1 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", TOTtime1[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nPARtime2 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", PARtime2[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nLTStime2 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", LTStime2[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nSUMtime2 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", SUMtime2[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nTOTtime2 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", TOTtime2[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */


    mexPrintf("\n\nTIME PERCENTAGE MATRICES\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nPARperc1 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", PARtime1[h][k]/TOTtime1[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nLTSperc1 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", LTStime1[h][k]/TOTtime1[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nSUMperc1 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", SUMtime1[h][k]/TOTtime1[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nPARperc2 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", PARtime2[h][k]/TOTtime2[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nLTSperc2 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", LTStime2[h][k]/TOTtime2[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    mexPrintf("\nSUMperc2 = [%c", ch);    for (k=0; k<NX; k++)    mexPrintf("%13d  ", NXval[k]);
    mexPrintf("= NXval\n");
    for (h=0; h<NT; h++)
    {   mexPrintf(  "            ");
        for (k=0; k<NX; k++)    mexPrintf("  %e", SUMtime2[h][k]/TOTtime2[h][k]);
        mexPrintf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    mexPrintf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */

}



void SEQ_PDE_skill_MAIN (double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol), /* user-defined function for LT samples */
                         unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol,
                         double *PARtime, double *LTStime, double *SUMtime, double *TOTtime)
{
    /* TALBOT SUITE DE FUNCTIONS */
    char *Talbot_fun[] = {"SEQ_TalbotSUM1_DE",  /* modified Talbot's method  */
                          "SEQ_TalbotSUM2_DE"}; /* classical Talbot's method */


    /* INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION */
    double sigma0;            /* abscissa of convergence                       */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    /*  Output value from function
            NsingsTOT: total number of singularities s_j
            unsigned int NsingsTOT = LTsings2(&Nsings, SINGS, MULT, &sigma0);
        is not used
    */
  /*unsigned int NsingsTOT = LTsings2(&Nsings, SINGS, MULT, &sigma0); */
    LTsings2(&Nsings, SINGS, MULT, &sigma0);


    /* INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM */
    double *t = linspace(NTval,Tmin,Tmax);

    /* INPUT DATA ABOUT ODE, COMING FROM APPLYING THE LAPLACE METHOD TO PDE */
    double *x = linspace(NXval,Xmin,Xmax);


    /** OUTPUT DATA **/
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot, *IFAIL; /* global and local error flags */
    NUMft = (double *)malloc(NXval*NTval*sizeof(double));
    if (NUMft == NULL)
         mexErrMsgTxt("\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n");

    IFAIL = (int *)malloc(NXval*NTval*sizeof(int));
    if (IFAIL == NULL)
         mexErrMsgTxt("\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n");


	/** VARIABLES OF THE main FUNCTION FOR TIMING **/
    double TOT_TIME;
    unsigned int NTIMES = 5, Nt;
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


    /* TALBOT SUITE DE: SKILL LEVEL */
    double CONLAM, CONSIG, CONNU; unsigned int NOPTS;
    double complex *S, *FS;
    IFAIL_tot = 0;
    int IFAIL_col = 0;
    unsigned int minNOPTS, maxNOPTS;

    switch (jFUN)
    {
    case 1 : /* modified Talbot's method for DE */

        *PARtime = 0.0;
        *LTStime = 0.0;
        *SUMtime = 0.0;
        for (Nt=0; Nt<NTIMES; Nt++) /* repeatitions */
        {
            /* 1.1) Compute Talbot's parameters at the middle point of the interval [Tmin, Tmax] */
            #ifdef __unix__
                clock_gettime(CLOCK_REALTIME, &ts1);
            #endif
            #if defined(_WIN64) || defined(_WIN32)
                QueryPerformanceCounter(&TICKS1);
            #endif
            /* ************************************************************************************** */
            COM_TalbotPAR (sigma0,(Tmin + Tmax)/2,tol,Nsings,SINGS,MULT, &CONLAM,&CONSIG,&CONNU,&NOPTS);
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
            *PARtime += TOT_TIME;

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
                mexErrMsgTxt("\n***   ERROR in SEQ_PDE_skill_MAIN: DYNAMIC ALLOCATION OF S IS FAILED. ***\n");

            double thetak,   piN = 4.0*atan(1.0)/NOPTS; /* pi over NOPTS (step) */
            unsigned int k;
            for (k=1; k<NOPTS; k++)
            {   thetak = k*piN;
                S[k] = CONSIG + CONLAM*thetak/tan(thetak) + I*CONLAM*CONNU*thetak; /* in math.h cot(x) doesn't exist */
            }
            S[0] = CONSIG + CONLAM;
            FS = (*LTsamples)(NXval,x,NOPTS,S,tol);
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
            *LTStime += TOT_TIME;

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
            *SUMtime += TOT_TIME;
        }

        *PARtime /= NTIMES; /* mean time for parameters */
        *LTStime /= NTIMES; /* mean time for LT samples */
        *SUMtime /= NTIMES; /* mean time for summations */
        *TOTtime = *PARtime + *LTStime + *SUMtime;
        break;


    case 2 : /* classical Talbot's method for DE */

        *PARtime = 0.0;
        *LTStime = 0.0;
        *SUMtime = 0.0;
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
                *PARtime += TOT_TIME;

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
                   mexErrMsgTxt("\n***   ERROR in SEQ_PDE_skill_MAIN: DYNAMIC ALLOCATION OF S IS FAILED. ***\n");

                double thetak,   piN = 4.0*atan(1.0)/NOPTS; /* pi over NOPTS (step) */
                unsigned int k;
                for (k=1; k<NOPTS; k++)
                {   thetak = k*piN;
                    S[k] = CONSIG + CONLAM*thetak/tan(thetak) + I*CONLAM*CONNU*thetak; /* in math.h cot(x) doesn't exist */
                }
                S[0] = CONSIG + CONLAM;
                FS = (*LTsamples)(NXval,x,NOPTS,S,tol);
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
                *LTStime += TOT_TIME;

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
                *SUMtime += TOT_TIME;
            } /* end-for j<NTval */
        } /* end for Nt<NTIMES */

        *PARtime /= NTIMES; /* mean time for parameters */
        *LTStime /= NTIMES; /* mean time for LT samples */
        *SUMtime /= NTIMES; /* mean time for summations */
        *TOTtime = *PARtime + *LTStime + *SUMtime;
        break;


    default :
            mexPrintf( "\n*** No function corresponds to the value %d ***", jFUN);
            exit(1);

    } /* END switch jFUN */

    /* PRINT RESULTS */
    mexPrintf("\nInverting the LT fun at NTval = %4d  values of t in [%5.2f, %5.2f]", NTval,Tmin,Tmax);
    mexPrintf("\n                 and at NXval = %4d  values of x in [%5.2f, %5.2f]", NXval,Xmin,Xmax);
    mexPrintf("\n                 input tolerance: tol = %5.2e\n", tol);
    mexPrintf("\n\t***   RESULTS OF SEQUENTIAL TALBOT SUITE DE [function = %s()]   ***\n", Talbot_fun[jFUN-1]);
    mexPrintf("\n\tExecution times:");
    mexPrintf("\n\tPAR time           = %e     for Talbot's parameters", *PARtime);
    mexPrintf("\n\tLT sample time     = %e     for LT samples of U(x,s) on Talbot's contour", *LTStime);
    mexPrintf("\n\tSUM time           = %e     for the Talbot-Clenshaw summation", *SUMtime);
    if (jFUN == 1)
        mexPrintf("\t[NOPTS = %d]", NOPTS);
    else
        mexPrintf("\t[%d <= NOPTS <= %d]", minNOPTS,maxNOPTS);
    mexPrintf("\n\tTotal elapsed time = %e\n", *TOTtime);

    free(NUMft); free(IFAIL);
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
