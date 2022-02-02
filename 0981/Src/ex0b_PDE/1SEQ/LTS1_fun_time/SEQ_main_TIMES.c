/********************        SEQ_main_TIMES.c        **********************
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *           DRIVER PROGRAM FOR THE APPLICATION TO PDE PROBLEMS           *
 *            OF THE SEQUENTIAL IMPLEMENTATION OF TALBOT SUITE            *
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *  >>>>>>>>>>>>>       VERSION 4.0   May 15th, 2016         <<<<<<<<<<<  *
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *                   AUTHOR: Mariarosaria Rizzardi                        *
 *                                                                        *
 *                   DiST   - "Parthenope" University, Naples (Italy)  *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 * REFERENCES                                                             *
 * ==========                                                             *
 * M. RIZZARDI: "Algorithm xxx: APPLICATION OF MODIFIED TALBOT'S METHOD   *
 *                              TO SOLVE DIFFERENTIAL PROBLEMS".          *
 *                    ACM Trans. Math. Softw., vol.##, no.#, month year,  *
 *                    pp. ##-##.                                          *
 *                                                                        *
 **************************************************************************
                               PDE EXAMPLE 0b
                          Example 2.2.6 (pag.88)
 Dean G. Duffy - Transform Methods for solving partial differential equations.
                 Chapman & Hall/CRC, 2004

 Application of Talbot's method to invert the following Laplace Transform

                U(eta,s) = exp(-eta*(1+sqrt(s+1.0)))*s/(s^2+9)^2

 U(X0,s) is a LT fun with poles of 2nd order at s=+/-3i and a branch point
 at s=-1, so that geometrical and accuracy parameters for Talbot's method
 have to be computed with the following informations on the singularities:
                    Nsings = 2;
                    Ssings = [-1 3i];
                    MULT   = [ 0  2];

**************************************************************************/

#include "../../../TalbotSuiteDE/COM/main_include.h"
#include "../../../TalbotSuiteDE/COM/COM_Talbot_pack.h"
#include "../../../TalbotSuiteDE/FUN/SEQ_Talbot_pack.h"       /* for Talbot Suite */
#include "../../../TalbotSuiteDE/FUN_DE/SEQ_Talbot_pack_DE.h" /* for Talbot Suite DE */

#include "../../problem_string.h"

#include <time.h>               /* for timing */

#if defined _WIN64 || defined _WIN32
    #include <windows.h>        /* for Windows high resolution timing */
#endif


/** global variable for eta in U(eta,s) [REQUIRED ONLY TO USE Talbot Suite] **/
double eta;


/** PROTOTYPES **/
unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0);
double *linspace (unsigned int Nv, double Vmin, double Vmax);


/* F(s): input parameter used by the Talbot Suite functions */
double complex LTfun(double complex s);
void SEQ_PDE_skill_MAIN (double complex (*LTpt)(double complex s), /* user-defined function for LT samples */
                         unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol, double *TOTtime);


/* U(x,s): input parameter used by the Talbot Suite DE functions */
double complex *SEQ_LTsamples_fun (unsigned int NXval, double X[], unsigned int NOPTS, double complex S[], double tol);
void SEQ_PDE_skill_MAIN_DE (double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol), /* user-defined function for LT samples */
                            unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol, double *TOTtime);




int main(int argc, char *argv[])
/*  To call:
                ./ex0.exe 1e-6

    Parameter to the main function:
        1) tol: tolerance on error
*/
{
    char *LTS = "LT samples computed by a function";

    double tol;                             /* tolerance            */
    if (argc > 1)
        tol = (double)atof(argv[1]);        /* 1st argument to main */
    else
        tol = 1e-8;                         /* default value        */


    /*   INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM              */
    const unsigned int NX=3, NT=3;
    unsigned int NTval[]={5,20,120},   NXval[]={5,20,120}; /* number of x, t values */
    double   Tmin=0.5,  Tmax=20.0,   Xmin=0.0,  Xmax=1.0;  /* endpoints */


    printf("%s", str);
    puts("\n\t   Ex. 0b: output from ./1SEQ/LTS1_fun/SEQ_main_TIMES.c");
    printf("\n\t%s\n", LTS);
    printf("\n\tTmin=%g;   Tmax=%g;   Xmin=%g;   Xmax=%g;   tol=%e\n", Tmin,Tmax,Xmin,Xmax,tol);
    puts("\n====================================================================================\n");


    /*  Mean execution time T(h,k) for NTval(h), NXval(k) */
    double TOTtime1[NT][NX],
           TOTtime2[NT][NX],
           TOTtime1_DE[NT][NX],
           TOTtime2_DE[NT][NX];


    unsigned int h, k; /* for-loop index */


    /* COMPUTE THE TIME T(h,k) */
    for (h=0; h<NT; h++)        /* Tval[h] */
        for (k=0; k<NX; k++)    /* Xval[k] */
        {
            SEQ_PDE_skill_MAIN (LTfun,1,NTval[h],Tmin,Tmax,NXval[k],Xmin,Xmax,tol,&TOTtime1[h][k]);
            SEQ_PDE_skill_MAIN (LTfun,1,NTval[h],Tmin,Tmax,NXval[k],Xmin,Xmax,tol,&TOTtime1[h][k]);
            /* chiama 2 volte Talbot1 per ammortizzare l'inizializzazione dei cicli for */
            SEQ_PDE_skill_MAIN (LTfun,2,NTval[h],Tmin,Tmax,NXval[k],Xmin,Xmax,tol,&TOTtime2[h][k]);

            SEQ_PDE_skill_MAIN_DE (SEQ_LTsamples_fun,1,NTval[h],Tmin,Tmax,NXval[k],Xmin,Xmax,tol,&TOTtime1_DE[h][k]);
            SEQ_PDE_skill_MAIN_DE (SEQ_LTsamples_fun,2,NTval[h],Tmin,Tmax,NXval[k],Xmin,Xmax,tol,&TOTtime2_DE[h][k]);

        }


    /* DISPLAY TIME'S MATRICES: 0 (LT function), 1 (LT ode.c), 2 (LT ode45.m) */
    char ch = '%';
    printf("\n\ttol = %e\n\n", tol);
/* -------------------------------------------------------------------------------------
        Talbot Suite
   ------------------------------------------------------------------------------------- */
    printf("\nTOTtime01 = [%c", ch);    for (k=0; k<NX; k++)    printf("%13d  ", NXval[k]);
    printf("= NXval\n");
    for (h=0; h<NT; h++)
    {   printf(  "            ");
        for (k=0; k<NX; k++)    printf("  %e", TOTtime1[h][k]);
        printf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    printf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    printf("\nTOTtime02 = [%c", ch);    for (k=0; k<NX; k++)    printf("%13d  ", NXval[k]);
    printf("= NXval\n");
    for (h=0; h<NT; h++)
    {   printf(  "            ");
        for (k=0; k<NX; k++)    printf("  %e", TOTtime2[h][k]);
        printf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    printf(  "            ];\n");
/* -------------------------------------------------------------------------------------
        Talbot Suite DE
   ------------------------------------------------------------------------------------- */
    printf("\nTOTtime01_DE = [%c", ch);    for (k=0; k<NX; k++)    printf("%13d  ", NXval[k]);
    printf("= NXval\n");
    for (h=0; h<NT; h++)
    {   printf(  "            ");
        for (k=0; k<NX; k++)    printf("  %e", TOTtime1_DE[h][k]);
        printf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    printf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */
    printf("\nTOTtime02_DE = [%c", ch);    for (k=0; k<NX; k++)    printf("%13d  ", NXval[k]);
    printf("= NXval\n");
    for (h=0; h<NT; h++)
    {   printf(  "            ");
        for (k=0; k<NX; k++)    printf("  %e", TOTtime2_DE[h][k]);
        printf("\t%c NTval = %3d\n", ch,NTval[h]);
    }
    printf(  "            ];\n");
/* ------------------------------------------------------------------------------------- */

    return 0;
}


void SEQ_PDE_skill_MAIN (double complex (*LTpt)(double complex s), /* user-defined function for LT samples */
                         unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol, double *TOTtime)
{
    /** TALBOT SUITE DE FUNCTIONS **/
    char *Talbot_fun[] = {"SEQ_Talbot1",  /* modified Talbot's method  */
                          "SEQ_Talbot2"}; /* classical Talbot's method */


    /** INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION **/
    double sigma0;            /* abscissa of convergence                       */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    /*  Output value from function LTsings2 is not used
            unsigned int NsingsTOT = LTsings2(&Nsings, SINGS, MULT, &sigma0);
            NsingsTOT: total number of singularities s_j
    */
    LTsings2(&Nsings, SINGS, MULT, &sigma0);


    /** INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM **/
    double *t = linspace(NTval,Tmin,Tmax);

    /** INPUT DATA ABOUT ODE, COMING FROM APPLYING THE LAPLACE METHOD TO PDE **/
    double *x = linspace(NXval,Xmin,Xmax);


    /** OUTPUT DATA **/
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot=0,   IFAIL_loc,   *IFAIL; /* global and local error flags */
    NUMft = (double *)malloc(NXval*NTval*sizeof(double));
    if (NUMft == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n"); exit(1);}

    IFAIL = (int *)malloc(NXval*NTval*sizeof(int));
    if (IFAIL == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n"); exit(1);}


	/** VARIABLES OF THE main FUNCTION FOR TIMING **/
    double TIME;
    unsigned int NTIMES = 5,  Nt, h;
    #ifdef __unix__
        printf("\nThis is Unix: %u bit C library\n", (unsigned int)sizeof(void*)*8 );
        struct timespec ts1, ts2; double T0, T1;
    #endif

    #if defined(_WIN32) || defined(_WIN64)
        printf("\nThis is Windows: %u bit C library\n", (unsigned int)sizeof(void*)*8);
        LARGE_INTEGER ticksPerSecond, TICKS1, TICKS2;
        QueryPerformanceFrequency(&ticksPerSecond); /* processor clock frequency */
        printf("\tCLOCK FREQUENCY by QueryPerformanceFrequency(): %g GHz\n\n",(double)ticksPerSecond.QuadPart/1e6);
    #endif



    /** TALBOT SUITE: user level **/

    switch (jFUN)
    {
    case 1 : /* MODIFIED TALBOT'S METHOD */

        *TOTtime = 0;
        for (Nt=0; Nt<NTIMES; Nt++) /* repeatitions */
        {
            #ifdef __unix__
                clock_gettime(CLOCK_REALTIME, &ts1);
            #endif
            #if defined(_WIN64) || defined(_WIN32)
                QueryPerformanceCounter(&TICKS1);
            #endif
            /* ************************************************************************************** */
            for (h=0; h<NXval; h++)
            {
                eta = x[h];
                IFAIL_loc = SEQ_Talbot1 ( LTpt, sigma0,NTval,t,tol,NUMft+h*NTval,IFAIL+h*NTval,Nsings,SINGS,MULT,Tmin,Tmax);
                IFAIL_tot = IFAIL_tot || IFAIL_loc;
            }
            /* ************************************************************************************** */
            #ifdef __unix__
                clock_gettime(CLOCK_REALTIME, &ts2);
                T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec      */
                T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
                TIME = T0 + T1*1.0e-9;
            #endif
            #if defined(_WIN64) || defined(_WIN32)
                QueryPerformanceCounter(&TICKS2);
                TIME = (double)(TICKS2.QuadPart - TICKS1.QuadPart)/(double)ticksPerSecond.QuadPart;
            #endif
            *TOTtime += TIME;
        }
        *TOTtime /= NTIMES; /* mean time */
        break;


    case 2 : /* CLASSICAL TALBOT'S METHOD */

        *TOTtime = 0;
        for (Nt=0; Nt<NTIMES; Nt++) /* repeatitions */
        {
            #ifdef __unix__
                clock_gettime(CLOCK_REALTIME, &ts1);
            #endif
            #if defined(_WIN64) || defined(_WIN32)
                QueryPerformanceCounter(&TICKS1);
            #endif
            /* ************************************************************************************** */
            for (h=0; h<NXval; h++)
            {
                eta = x[h];
                IFAIL_loc = SEQ_Talbot2 ( LTpt, sigma0,NTval,t,tol,NUMft+h*NTval,IFAIL+h*NTval,Nsings,SINGS,MULT);
                IFAIL_tot = IFAIL_tot || IFAIL_loc;
            }
            /* ************************************************************************************** */
            #ifdef __unix__
                clock_gettime(CLOCK_REALTIME, &ts2);
                T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec      */
                T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
                TIME = T0 + T1*1.0e-9;
            #endif
            #if defined(_WIN64) || defined(_WIN32)
                QueryPerformanceCounter(&TICKS2);
                TIME = (double)(TICKS2.QuadPart - TICKS1.QuadPart)/(double)ticksPerSecond.QuadPart;
            #endif
            *TOTtime += TIME;
        } /* end for Nt<NTIMES */
        *TOTtime /= NTIMES; /* mean time */
        break;


    default :
        fprintf(stderr, "\n*** No function corresponds to the value %d ***", jFUN);
        exit(1);

    } /* END switch jFUN */

    /* PRINT RESULTS */
    printf("\nInverting the LT fun at NTval = %4d  values of t in [%5.2f, %5.2f]", NTval,Tmin,Tmax);
    printf("\n                 and at NXval = %4d  values of x in [%5.2f, %5.2f]", NXval,Xmin,Xmax);
    printf("\n                 input tolerance: tol = %5.2e\n", tol);
    printf("\n\t***   RESULTS OF SEQUENTIAL TALBOT SUITE [function = %s()]   ***\n", Talbot_fun[jFUN-1]);
    printf("\n\tTotal elapsed time = %e\n", *TOTtime);

    free(NUMft); free(IFAIL);
}



void SEQ_PDE_skill_MAIN_DE (double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol), /* user-defined function for LT samples */
                            unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol, double *TOTtime)
{
    /** TALBOT SUITE DE FUNCTIONS **/
    char *Talbot_fun[] = {"SEQ_TalbotSUM1_DE",  /* modified Talbot's method  */
                          "SEQ_TalbotSUM2_DE"}; /* classical Talbot's method */


    /** INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION **/
    double sigma0;            /* abscissa of convergence                       */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    /*  Output value from function LTsings2 is not used
            unsigned int NsingsTOT = LTsings2(&Nsings, SINGS, MULT, &sigma0);
            NsingsTOT: total number of singularities s_j
    */
    LTsings2(&Nsings, SINGS, MULT, &sigma0);


    /** INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM **/
    double *t = linspace(NTval,Tmin,Tmax);

    /** INPUT DATA ABOUT ODE, COMING FROM APPLYING THE LAPLACE METHOD TO PDE **/
    double *x = linspace(NXval,Xmin,Xmax);


    /** OUTPUT DATA **/
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot,  *IFAIL; /* global and local error flags */
    NUMft = (double *)malloc(NXval*NTval*sizeof(double));
    if (NUMft == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n"); exit(1);}

    IFAIL = (int *)malloc(NXval*NTval*sizeof(int));
    if (IFAIL == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n"); exit(1);}


	/** VARIABLES OF THE main FUNCTION FOR TIMING **/
    double TIME;
    unsigned int NTIMES = 5, Nt;
    #ifdef __unix__
        printf("\nThis is Unix: %u bit C library\n", (unsigned int)sizeof(void*)*8 );
        struct timespec ts1, ts2; double T0, T1;
    #endif

    #if defined(_WIN32) || defined(_WIN64)
        printf("\nThis is Windows: %u bit C library\n", (unsigned int)sizeof(void*)*8);
        LARGE_INTEGER ticksPerSecond, TICKS1, TICKS2;
        QueryPerformanceFrequency(&ticksPerSecond); /* processor clock frequency */
        printf("\tCLOCK FREQUENCY by QueryPerformanceFrequency(): %g GHz\n\n",(double)ticksPerSecond.QuadPart/1e6);
    #endif


    /** TALBOT SUITE DE: user level **/

    switch (jFUN)
    {
    case 1 : /* modified Talbot's method for DE */

        *TOTtime = 0.0;
        for (Nt=0; Nt<NTIMES; Nt++) /* repeatitions */
        {

            #ifdef __unix__
                clock_gettime(CLOCK_REALTIME, &ts1);
            #endif
            #if defined(_WIN64) || defined(_WIN32)
                QueryPerformanceCounter(&TICKS1);
            #endif
            /* ************************************************************************************** */
            IFAIL_tot = SEQ_Talbot1_DE (LTsamples, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax);
            /* ************************************************************************************** */
            #ifdef __unix__
                clock_gettime(CLOCK_REALTIME, &ts2);
                T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec      */
                T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
                TIME = T0 + T1*1.0e-9;
            #endif
            #if defined(_WIN64) || defined(_WIN32)
                QueryPerformanceCounter(&TICKS2);
                TIME = (double)(TICKS2.QuadPart - TICKS1.QuadPart)/(double)ticksPerSecond.QuadPart;
            #endif
            *TOTtime += TIME;
        }
        *TOTtime /= NTIMES; /* mean time */
        break;


    case 2 : /* classical Talbot's method for DE */

        *TOTtime = 0.0;
        for (Nt=0; Nt<NTIMES; Nt++) /* repeatitions */
        {
            #ifdef __unix__
                clock_gettime(CLOCK_REALTIME, &ts1);
            #endif
            #if defined(_WIN64) || defined(_WIN32)
                QueryPerformanceCounter(&TICKS1);
            #endif
            /* ************************************************************************************** */
            IFAIL_tot = SEQ_Talbot2_DE (LTsamples, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT);
            /* ************************************************************************************** */
            #ifdef __unix__
                clock_gettime(CLOCK_REALTIME, &ts2);
                T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec      */
                T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
                TIME = T0 + T1*1.0e-9;
            #endif
            #if defined(_WIN64) || defined(_WIN32)
                QueryPerformanceCounter(&TICKS2);
                TIME = (double)(TICKS2.QuadPart - TICKS1.QuadPart)/(double)ticksPerSecond.QuadPart;
            #endif
            *TOTtime += TIME;
        }
        *TOTtime /= NTIMES; /* mean time */
        break;


    default :
        fprintf(stderr, "\n*** No function corresponds to the value %d ***", jFUN);
        exit(1);

    } /* END switch jFUN */

    /* PRINT RESULTS */
    printf("\nInverting the LT fun at NTval = %4d  values of t in [%5.2f, %5.2f]", NTval,Tmin,Tmax);
    printf("\n                 and at NXval = %4d  values of x in [%5.2f, %5.2f]", NXval,Xmin,Xmax);
    printf("\n                 input tolerance: tol = %5.2e\n", tol);
    printf("\n\t***   RESULTS OF SEQUENTIAL TALBOT SUITE DE [function = %s()]   ***\n", Talbot_fun[jFUN-1]);
    printf("\n\tTotal elapsed time = %e\n", *TOTtime);

    free(NUMft); free(IFAIL);
}



double *linspace(unsigned int Nv, double Vmin, double Vmax)
/** Compute a dynamic double array of equispace values **/
{
    double *v = (double *)malloc(Nv*sizeof(double));
    if (v == NULL)
        {fprintf(stderr, "\n***   ERROR IN linspace(): DYNAMIC ALLOCATION OF THE ARRAY IS FAILED. ***\n"); exit(1);}

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
