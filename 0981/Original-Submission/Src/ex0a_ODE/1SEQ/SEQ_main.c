/************************        SEQ_main.c       *************************
 *                                                                        *
 *                                                                        *
 *                        SEQUENTIAL TALBOT SUITE DE                      *
 *                                                                        *
 *                       DRIVER PROGRAM FOR EXAMPLE 0a                    *
 *                                                                        *
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
 * M. RIZZARDI: "Algorithm xxx: APPLICATION OF MODIFIED TALBOT'S METHOD   *
 *                              TO SOLVE DIFFERENTIAL PROBLEMS".          *
 *                                 ACM TRANS. MATH. SOFTWARE, VOL. xx,    *
 *                                 NO. x, month year, (PP ##).            *
 *                                                                        *
 **************************************************************************/

#include "../../TalbotSuiteDE/COM/main_include.h"          /* for shared #include                             */
#include "../../TalbotSuiteDE/COM/COM_Talbot_pack.h"       /* for Talbot Suite's shared functions             */
#include "../../TalbotSuiteDE/COM_DE/COM_Talbot_pack_DE.h" /* for correction to NOPTS                         */
#include "./SEQ_Talbot_pack.h" /* for Talbot Suite's sequential implementation with no N correction */
#include "../../TalbotSuiteDE/FUN_DE/SEQ_Talbot_pack_DE.h" /* for Talbot Suite DE's sequential implementation */
#include "../problem_string.h"

#include <time.h>                   /* for timing               */

#if defined _WIN64 || defined _WIN32
    #include <windows.h>            /* for timing under Windows */
#endif


/*  Inverse LT function: not required
    =================================
    This function has been included to compare the numerical approximations
    to the corresponding true results
 */
double ILTfun (double t);


/* Abscissa of convergence and singularities of F(s) with their polar multiplicities: required */
unsigned int LTsings (unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0);
double complex LTfun (double complex s);
double complex *SEQ_LTsamples_fun (unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol);


int main(int argc, char *argv[])
/*  ./ex0.exe test

    test: 1st argument to main function
         switch among 3 test examples
           test = 1      20 t in [1000, 3000]
                = 2     120 t in [1000, 3000]
                = 3     120 t in [ 100,  500]
                = 4      20 t in [ 100,  500]
 */
{
    unsigned int test = 0;
    if (argc > 1)
        test = (unsigned int)atoi(argv[1]); /* 1st argument to main */
    else
        test = 1;                           /* default value        */



    char *LTS[] = {"LTfun",              /* modified  method Talbot Suite    */
                   "LTfun",              /* classical method Talbot Suite    */
                   "SEQ_LTsamples_fun",  /* modified  method Talbot Suite DE */
                   "SEQ_LTsamples_fun"}; /* classical method Talbot Suite DE */

    char *Talbot_fun[] = {"SEQ_Talbot1",     /* modified  method Talbot Suite    */
                          "SEQ_Talbot2",     /* classical method Talbot Suite    */
                          "SEQ_Talbot1_DE",  /* modified  method Talbot Suite DE */
                          "SEQ_Talbot2_DE"}; /* classical method Talbot Suite DE */


    /* INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION */
    double sigma0;            /* abscissa of convergence                       */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    LTsings(&Nsings, SINGS, MULT, &sigma0);


    unsigned int k; /* for-loop index */
    unsigned int PRINTflag = 1; /* false: do not print accuracy */


    /*   INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM        */
    double tol = 1e-12;    /* tolerance                       */
    unsigned int NXval=1;  /* dummy number of x values        */
    double *x;             /* pointer to dummy x variable     */
    unsigned int NTval;    /* number of t values              */
    double Tmin, Tmax;     /* endpoints of the interval for t */
    double *t;             /* t array                         */

    switch (test)
    {
        case 1 : /*  20 t in [1000, 3000] */
            NTval=20;  Tmin=1000;  Tmax=3000;
            break;
        case 2 : /* 120 t in [1000, 3000] */
            NTval=120;  Tmin=1000;  Tmax=3000;
            break;
        case 3 : /* 120 t in [ 100,  500] */
            NTval=120;  Tmin=100;  Tmax=500;
            break;
        case 4 : /*  20 t in [ 100,  500] */
            NTval=20;  Tmin=100;  Tmax=500;
            break;
        default :
            fprintf(stderr, "\n***   ERROR: The parameter to main() must be 1, 2 or 3. ***\n");
            exit(1);
    }

    printf("%s", str);
    puts("\n\t   Ex. 0a: output from ./1SEQ/LTS1_fun/SEQ_main.c");
    printf("\n\t%d t in [%g, %g],    tol=%e", NTval,Tmin,Tmax,tol);
    puts("\n====================================================================================\n");


    /*   OUTPUT DATA                                         */
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot, *IFAIL; /* global and local error flags */


    /* DUMMY x VALUE */
    x = (double *)malloc(NXval*sizeof(double));
    if (x == NULL)  { fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF x IS FAILED. ***\n"); exit(1); }
    *x = 0;


    /* COMPUTE THE t VALUES */
    t = (double *)malloc(NTval*sizeof(double));
    if (t == NULL)  { fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF t IS FAILED. ***\n"); exit(1); }
    double t_step;
    if (NTval>1)
        /* several values */
        t_step = (Tmax-Tmin)/(NTval-1);
    else
        /* a single value */
        t_step = 0.0;

    for (k=0; k<NTval; k++)
        *(t+k) = Tmin + k*t_step;


    /* ALLOCATE ALGORITHM'S OUTPUT ARRAYS */
    NUMft = (double *)malloc(NTval*sizeof(double));
    if (NUMft == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n"); exit(1);}

    IFAIL = (int *)malloc(NTval*sizeof(int));
    if (IFAIL == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n"); exit(1);}


    /* VARIABLES OF THE main FUNCTION */
    double ft, RELERR, ABSERR; char ALFA;


    /* FOR TIMING [both Windows and Linux] */
    double TOT_TIME, MEAN_TIME;
    unsigned int NTIMES = 5, Nt;
    #ifdef __unix__
        printf("This is Unix: %u bit C library\n", (unsigned int)sizeof(void*)*8 );
        struct timespec ts1, ts2; double T0, T1;
    #endif // __unix__

    #if defined(_WIN32) || defined(_WIN64)
        printf("This is Windows: %u bit C library\n", (unsigned int)sizeof(void*)*8);
        LARGE_INTEGER ticksPerSecond, TICKS1, TICKS2;
        QueryPerformanceFrequency(&ticksPerSecond); // processor clock frequency
        printf("\tCLOCK FREQUENCY by QueryPerformanceFrequency(): %g GHz\n",(double)ticksPerSecond.QuadPart/1e6);
    #endif


    printf("\n*******   Inverting the LT fun F(s) to approximate f(t) at NTval = %d values of t in [%5.2f, %5.2f]   *******\n\n", NTval,Tmin,Tmax);

    /* SELECTION OF FUNCTION TO BE TESTED */
    unsigned int jFUN; /* switch between SEQ_Talbot1    (1),
                                         SEQ_Talbot2    (2),
                                         SEQ_Talbot1_DE (3),
                                         SEQ_Talbot2_DE (4) */
    for (jFUN=1; jFUN<=4; jFUN++)
    {
        /* CALL THE TALBOT SUITE FUNCTION  */
        switch (jFUN)
        {
        case 1 : /* (modified method of Talbot Suite)  */
            MEAN_TIME = 0;
            for (Nt=0; Nt<NTIMES; Nt++)
            {
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts1);
                #endif // __unix__
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS1);
                #endif
                        IFAIL_tot = SEQ_Talbot1 ( LTfun,sigma0,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax );
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts2);
                    T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec      */
                    T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
                    TOT_TIME = T0 + T1*1.0e-9;
                #endif // __unix__
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS2);
                    TOT_TIME = (double)(TICKS2.QuadPart - TICKS1.QuadPart)/(double)ticksPerSecond.QuadPart;
                #endif
                MEAN_TIME += TOT_TIME;
            }
            MEAN_TIME /= NTIMES;
            break;

        case 2 : /*  (classical method of Talbot Suite) */
            MEAN_TIME = 0;
            for (Nt=0; Nt<NTIMES; Nt++)
            {
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts1);
                #endif // __unix__
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS1);
                #endif
                    IFAIL_tot = SEQ_Talbot2 ( LTfun,sigma0,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT );
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts2);
                    T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec      */
                    T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
                    TOT_TIME = T0 + T1*1.0e-9;
                #endif // __unix__
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS2);
                    TOT_TIME = (double)(TICKS2.QuadPart - TICKS1.QuadPart)/(double)ticksPerSecond.QuadPart;
                #endif
                MEAN_TIME += TOT_TIME;
            }
            MEAN_TIME /= NTIMES;
            break;

        case 3 : /* (modified method of Talbot Suite DE)  */
            MEAN_TIME = 0;
            for (Nt=0; Nt<NTIMES; Nt++)
            {
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts1);
                #endif // __unix__
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS1);
                #endif
                        IFAIL_tot = SEQ_Talbot1_DE ( SEQ_LTsamples_fun, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax);
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts2);
                    T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec      */
                    T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
                    TOT_TIME = T0 + T1*1.0e-9;
                #endif // __unix__
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS2);
                    TOT_TIME = (double)(TICKS2.QuadPart - TICKS1.QuadPart)/(double)ticksPerSecond.QuadPart;
                #endif
                MEAN_TIME += TOT_TIME;
            }
            MEAN_TIME /= NTIMES;
            break;

        case 4 : /*  (classical method of Talbot Suite DE) */
            MEAN_TIME = 0;
            for (Nt=0; Nt<NTIMES; Nt++)
            {
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts1);
                #endif // __unix__
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS1);
                #endif
                    IFAIL_tot = SEQ_Talbot2_DE ( SEQ_LTsamples_fun, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT);
                #ifdef __unix__
                    clock_gettime(CLOCK_REALTIME, &ts2);
                    T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec      */
                    T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
                    TOT_TIME = T0 + T1*1.0e-9;
                #endif // __unix__
                #if defined(_WIN64) || defined(_WIN32)
                    QueryPerformanceCounter(&TICKS2);
                    TOT_TIME = (double)(TICKS2.QuadPart - TICKS1.QuadPart)/(double)ticksPerSecond.QuadPart;
                #endif
                MEAN_TIME += TOT_TIME;
            }
            MEAN_TIME /= NTIMES;
            break;

        default :
            fprintf(stderr, "\n*** No function corresponds to the value %d ***", jFUN);
            exit(1);

        } /* END switch jFUN */

        /* PRINT RESULTS */
        (jFUN == 1  ||  jFUN == 2) ? printf("\n\t***   RESULTS OF SEQUENTIAL TALBOT SUITE function:      %s()   ***", Talbot_fun[jFUN-1])
                                   : printf("\n\t***   RESULTS OF SEQUENTIAL TALBOT SUITE DE function:   %s()   ***", Talbot_fun[jFUN-1]);
        printf("\n\t      user defined function to compute the LT samples:  %s()", LTS[jFUN-1]);
        printf("\n\t      input tolerance: tol = %5.2e\t\tmean elapsed time = %e", tol,MEAN_TIME);

        if ( jFUN%2 ) /* jFUN = 1, 3: only for the modified method */
        {
            double CONLAM, CONSIG, CONNU; unsigned int NOPTS;
            COM_TalbotPAR (sigma0,(Tmin+Tmax)/2,tol,Nsings,SINGS,MULT, &CONLAM,&CONSIG,&CONNU,&NOPTS);
            if (jFUN == 3) /* Talbot Suite DE */
                NOPTS = COM_TalbotNcorr(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4);
            printf("\n\t      NOPTS = %u\n", NOPTS);
        }
        else
            printf("\n");

        if (PRINTflag  &&  NTval < 50)
        {
            if (!IFAIL_tot)
            {
                printf("\n       T       F EXACT          F APPROX        ABS ERR         REL ERR      TYPE    IFAIL_tot = %d (no local error)", IFAIL_tot);
                /* for each t, compute f(t) and errors, then write results */
                for (k=0; k<NTval; k++)
                {
                    ft = ILTfun( t[k] );
                    if (fabs(ft)> 1.0)
                        ALFA='R';
                    else
                        ALFA='A';
                    ABSERR = fabs( ft - NUMft[k] );
                    if (fabs(ft)<LDBL_MIN) /* underflow */
                        printf("\n %10.2f   %+e   %+e   %13.4e         -          %c", t[k],ft,NUMft[k],ABSERR,ALFA);
                    else
                    {
                        RELERR = ABSERR/fabs(ft);
                        printf("\n %10.2f   %+e   %+e   %e   %e    %c", t[k],ft,NUMft[k],ABSERR,RELERR,ALFA);
                    }
                }
            }
            else
            {
                puts("\n       T       F EXACT          F APPROX        ABS ERR         REL ERR      TYPE    IFAIL");
                /* for each t, compute f(t) and errors, then write results and error flags */
                for (k=0; k<NTval; k++)
                {
                    ft = ILTfun( t[k] );
                    if (fabs(ft)> 1.0)
                        ALFA='R';
                    else
                        ALFA='A';
                    ABSERR = fabs( ft - NUMft[k] );
                    if (fabs(ft)<LDBL_MIN) /* underflow */
                        printf("\n %10.2f   %+e   %+e   %13.4e         -          %c   %+3d", t[k],ft,NUMft[k],ABSERR,ALFA,IFAIL[k]);
                    else
                    {
                        RELERR = ABSERR/fabs(ft);
                        printf("\n %10.2f   %+e   %+e   %e   %e    %c   %+3d", t[k],ft,NUMft[k],ABSERR,RELERR,ALFA,IFAIL[k]);
                    }
                }
            }
        }
        puts("\n\n");
    } /* END for (jFUN) */

    return 0;
}
