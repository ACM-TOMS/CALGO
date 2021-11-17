/************************        OMP_main.c       *************************
 *                                                                        *
 *                                                                        *
 *             OpenMP-based parallel version of TALBOT SUITE DE           *
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
 * M. RIZZARDI: "Algorithm xxx: TALBOT SUITE DE: APPLICATION OF MODIFIED  *
 *                              TALBOT'S METHOD TO SOLVE DIFFERENTIAL     *
 *                              PROBLEMS".                                *
 *                                 ACM TRANS. MATH. SOFTWARE, VOL. xx,    *
 *                                 NO. x, month year, (PP ##).            *
 *                                                                        *
 **************************************************************************/

#include "../../TalbotSuiteDE/COM/main_include.h"          /* for shared #include                             */
#include "../../TalbotSuiteDE/COM/COM_Talbot_pack.h"       /* for Talbot Suite's shared functions             */
#include "../../TalbotSuiteDE/COM_DE/COM_Talbot_pack_DE.h" /* for correction to NOPTS                         */
#include "../../TalbotSuiteDE/FUN_DE/OMP_Talbot_pack_DE.h" /* for Talbot Suite's OMP parallel implementation  */
#include "../../TalbotSuiteDE/FUN_DE/SEQ_Talbot_pack_DE.h" /* for Talbot Suite's  sequential  implementation  */
#include "../problem_string.h"

#include <time.h>                   /* for timing               */

#if defined _WIN64 || defined _WIN32
    #include <windows.h>            /* for timing under Windows */
#endif

#ifdef _OPENMP
    #include <omp.h>
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
double complex *OMP_LTsamples_fun (unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS);

void display_errors(unsigned int PRINTflag, unsigned int NTval, int IFAIL_tot, double t[], double NUMft[], int IFAIL[]);
double *linspace (unsigned int Nv, double Vmin, double Vmax);


int main(int argc, char *argv[])
/* ./ex0.exe test jFUN

    test: 1st argument to main function
          switch among 3 test examples
               test = 1      20 t in [1000, 3000]
                    = 2     120 t in [1000, 3000]
                    = 3     120 t in [ 100,  500]

    jFUN: 2nd argument to main function
          switch among Talbot Suite DE functions
               jFUN = 0 : SEQ_Talbot1_DE  - sequential modified Talbot's method
                    = 1 : OMP_Talbot11_DE - OpenMP coarse-grained parallel modified Talbot's method
                    = 2 : OMP_Talbot12_DE - OpenMP fine-grained parallel modified Talbot's method
                    = 3 : OMP_Talbot13_DE - OpenMP nested coarse/fine-grained parallel modified Talbot's method
                    = 4 : SEQ_Talbot2_DE  - sequential classical Talbot's method
*/
{
    unsigned int test, jFUN;
    if (argc > 2)
    {   test = (unsigned int)atoi(argv[1]); /* 1st argument to main */
        jFUN = (unsigned int)atoi(argv[2]); /* 2nd argument to main */
    }
    else
    {   test = 3;                           /* default values       */
        jFUN = 1;
    }

    char *LTS = "OMP_LTsamples_fun";  /* OpenMP parallel LT samples */

    char *Talbot_fun[] = {"SEQ_Talbot1_DE",   /* 0: sequential modified Talbot's method */
                          "OMP_Talbot11_DE",  /* 1: OpenMP coarse grained parallel modified Talbot's method */
                          "OMP_Talbot12_DE",  /* 2: OpenMP fine grained parallel modified Talbot's method */
                          "OMP_Talbot13_DE",  /* 3: OpenMP nested coarse/fine grained parallel modified Talbot's method */
                          "SEQ_Talbot2_DE"}; /* 4: sequential classical Talbot's method */


    unsigned int k; /* for-loop index */
    unsigned int PRINTflag = 1; /* false: do not print accuracy
                                   true : print absolute error */


    /* INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION */
    double sigma0;            /* abscissa of convergence                       */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    LTsings(&Nsings, SINGS, MULT, &sigma0);


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
        default :
            fprintf(stderr, "\n***   ERROR: The 1st parameter to main() must be 1, 2 or 3. ***\n");
            exit(1);
    }

    printf("%s", str);
    puts("\n\t   Ex. 0a: output from ./2PAR/OMP_main.c");
    printf("\t%s\n", LTS);
    printf("\n\t%d t in [%g, %g],    tol=%e", NTval,Tmin,Tmax,tol);
    puts("\n====================================================================================\n");


    /*   OUTPUT DATA                                       */
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot, *IFAIL; /* global and local error flags */


    /* DUMMY x VALUE */
    x = (double *)malloc(NXval*sizeof(double));
    if (x == NULL)  { fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF x IS FAILED. ***\n"); exit(1); }
    *x = 0;


    /* COMPUTE THE t VALUES */
    t = linspace(NTval,Tmin,Tmax);


    /* ALLOCATE ALGORITHM'S OUTPUT ARRAYS */
    NUMft = (double *)malloc(NTval*sizeof(double));
    if (NUMft == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n"); exit(1);}

    IFAIL = (int *)malloc(NTval*sizeof(int));
    if (IFAIL == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n"); exit(1);}


    /* FOR TIMING [both Windows and Linux] */
    double TOT_TIME, MEAN_TIME;
    unsigned int NTIMES = 5, Nt;
    #ifdef __unix__
        printf("This is Unix: %u bit C library\n", (unsigned int)sizeof(void*)*8 );
        struct timespec ts1, ts2; double T0, T1;
    #endif

    #if defined(_WIN32) || defined(_WIN64)
        printf("This is Windows: %u bit C library\n", (unsigned int)sizeof(void*)*8);
        LARGE_INTEGER ticksPerSecond, TICKS1, TICKS2;
        QueryPerformanceFrequency(&ticksPerSecond); /* processor clock frequency */
        printf("\tCLOCK FREQUENCY by QueryPerformanceFrequency(): %g GHz\n",(double)ticksPerSecond.QuadPart/1e6);
    #endif


    printf("\n*******   Inverting the LT fun F(s) to approximate f(t) at NTval = %d values of t in [%5.2f, %5.2f]   *******", NTval,Tmin,Tmax);
    printf("\n          user defined function to compute the LT samples:  %s(),   input tolerance: tol = %5.2e\n\n", LTS,tol);


    int THREADS, maxTHREADS;
    #ifdef _OPENMP
        printf("\n\t[OpenMP 4.0] omp_get_num_procs() = %d,   omp_get_max_threads() = %d,   omp_get_num_threads() = %d\n", omp_get_num_procs(),omp_get_max_threads(),omp_get_num_threads() );
        maxTHREADS = omp_get_num_procs(); /* OpenMP 4.0 */
    #else
        printf("OpenMP not defined!\n");
        maxTHREADS = 1;
    #endif

    char ch='%';
    printf("\n\t***   RESULTS OF TALBOT SUITE DE function:   %s()   ***", Talbot_fun[jFUN]);
    if (jFUN == 1 || jFUN == 2 || jFUN == 3)
    {
        double CONLAM, CONSIG, CONNU; unsigned int NOPTS;
        COM_TalbotPAR (sigma0,(Tmin+Tmax)/2,tol,Nsings,SINGS,MULT, &CONLAM,&CONSIG,&CONNU,&NOPTS);
        NOPTS = COM_TalbotNcorr(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4);
        printf("\n\t      NOPTS = %u", NOPTS);
    }



    if (jFUN == 3)
    {   /* OpenMP nested coarse/fine grained parallel modified Talbot's method */
        int thrds1, thrds2, sqrMAX;
        #ifdef _OPENMP
            printf("\n\n\tCHECK: NESTED PARALLELISM is %s\n", omp_get_nested() ? "supported" : "not supported");
            if ( !omp_get_nested() )
            {   omp_set_nested(1); /* true */
                printf("\t                          Now it is %s\n", omp_get_nested() ? "supported" : "not supported");
            }
        #endif
        sqrMAX = ceil(sqrt((double)maxTHREADS)); /* for nested parallelism */
        printf("\nMEAN ELAPSED TIME:\n");
        printf("%c%15d", ch,1);
        for (k=2; k<=sqrMAX+2; k++)    printf("%16d", k);
        printf("  =    number of inner threads ( fine-grain  parallelism)\n");

        for (thrds1=1;  thrds1<=sqrMAX+2;  thrds1++) /* number of outer parallel threads */
        {
            for (thrds2=1;  thrds2<=sqrMAX+2;  thrds2++) /* number of inner parallel threads */
            {
                MEAN_TIME = NAN;
                if ( thrds1*thrds2 > maxTHREADS  &&  ( thrds1 > sqrMAX  ||  thrds2 > sqrMAX ) )
                {
                    printf("             nan");
                    continue;
                }
                else
                {
                    /* SET THE NUMBER OF OpenMP PARALLEL THREADS */
                    #ifdef _OPENMP
                        omp_set_dynamic(0); /* static scheduling */
                        omp_set_num_threads(thrds1*thrds2);
                    #endif
                    MEAN_TIME = 0;
                    for (Nt=0; Nt<NTIMES; Nt++)
                    {
                        #ifdef __unix__
                            clock_gettime(CLOCK_REALTIME, &ts1);
                        #endif
                        #if defined(_WIN64) || defined(_WIN32)
                            QueryPerformanceCounter(&TICKS1);
                        #endif
                                IFAIL_tot = OMP_Talbot13_DE ( OMP_LTsamples_fun, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax,thrds1,thrds2);
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
                }
                printf("%16.4e", MEAN_TIME); /* PRINT RESULTS */
            }
            printf("  %c %2d%s", ch,thrds1, thrds1 == 1 ? " number of outer threads (coarse-grain parallelism)\n" : "\n");
        }
        puts("");
    }
    else
    {
        printf("\n\tMEAN ELAPSED TIME:\n");
        /* Loop on number of threads */
        for ( THREADS=1;  THREADS<=maxTHREADS;  THREADS++ )
        {
            /* SET THE NUMBER OF OpenMP PARALLEL THREADS */
            #ifdef _OPENMP
                omp_set_dynamic(0); /* static scheduling */
                omp_set_num_threads(THREADS);
            #endif

            /* CALL THE TALBOT SUITE FUNCTION  */
            switch (jFUN) {

            case 0 : /* sequential modified Talbot's method */
                MEAN_TIME = 0;
                for (Nt=0; Nt<NTIMES; Nt++)
                {
                    #ifdef __unix__
                        clock_gettime(CLOCK_REALTIME, &ts1);
                    #endif
                    #if defined(_WIN64) || defined(_WIN32)
                        QueryPerformanceCounter(&TICKS1);
                    #endif
                            IFAIL_tot = SEQ_Talbot1_DE ( SEQ_LTsamples_fun, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax);
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
                break;

            case 1 : /* OpenMP coarse grained parallel modified Talbot's method */
                MEAN_TIME = 0;
                for (Nt=0; Nt<NTIMES; Nt++)
                {
                    #ifdef __unix__
                        clock_gettime(CLOCK_REALTIME, &ts1);
                    #endif
                    #if defined(_WIN64) || defined(_WIN32)
                        QueryPerformanceCounter(&TICKS1);
                    #endif
                            IFAIL_tot = OMP_Talbot11_DE ( OMP_LTsamples_fun, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax,THREADS);
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
                break;

            case 2 : /* OpenMP fine grained parallel modified Talbot's method */
                MEAN_TIME = 0;
                for (Nt=0; Nt<NTIMES; Nt++)
                {
                    #ifdef __unix__
                        clock_gettime(CLOCK_REALTIME, &ts1);
                    #endif
                    #if defined(_WIN64) || defined(_WIN32)
                        QueryPerformanceCounter(&TICKS1);
                    #endif
                            IFAIL_tot = OMP_Talbot12_DE ( OMP_LTsamples_fun, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax,THREADS);
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
                break;

            case 4 : /*  sequential classical Talbot's method */
                MEAN_TIME = 0;
                for (Nt=0; Nt<NTIMES; Nt++)
                {
                    #ifdef __unix__
                        clock_gettime(CLOCK_REALTIME, &ts1);
                    #endif
                    #if defined(_WIN64) || defined(_WIN32)
                        QueryPerformanceCounter(&TICKS1);
                    #endif
                        IFAIL_tot = SEQ_Talbot2_DE ( SEQ_LTsamples_fun, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT);
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
                break;

            default :
                fprintf(stderr, "\n*** No function corresponds to the value %d ***", jFUN);
                exit(1);

            } /* END switch jFUN */

            /* PRINT RESULTS */
            printf("\t%e\t%c number of threads = %2d\n", MEAN_TIME,ch,THREADS);

        } /* END for (THREADS) */
    }

    display_errors(PRINTflag,NTval,IFAIL_tot,t,NUMft,IFAIL);
    puts("\n\n");
    return 0;
}


void display_errors(unsigned int PRINTflag, unsigned int NTval, int IFAIL_tot, double t[], double NUMft[], int IFAIL[])
{
    char ALFA;
    unsigned int k;
    double ft, RELERR, ABSERR;

    if (PRINTflag  &&  NTval < 50)
    {
        if (!IFAIL_tot)
        {
            printf("\n       T       F EXACT          F APPROX        ABS ERR         REL ERR      ERR TYPE    IFAIL_tot = %d (no local error)", IFAIL_tot);
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
    puts("\n\n");
    }
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
