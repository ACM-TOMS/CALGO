/*********************        OMP_main_TIMES.c       **********************
 *                                                                        *
 *                                                                        *
 *                  OpenMP-BASED PARALLEL TALBOT SUITE DE                 *
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

#include "../../../TalbotSuiteDE/COM/main_include.h"          /* for shared #include                             */
#include "../../../TalbotSuiteDE/COM/COM_Talbot_pack.h"       /* for Talbot Suite's shared functions             */
#include "../../../TalbotSuiteDE/COM_DE/COM_Talbot_pack_DE.h" /* for correction to NOPTS                         */
#include "../../../TalbotSuiteDE/FUN_DE/OMP_Talbot_pack_DE.h" /* for Talbot Suite's OMP parallel implementation  */
#include "../../problem_string.h"

#include <time.h>                   /* for timing                               */

#if defined _WIN64 || defined _WIN32
    #include <windows.h>        /* for LARGE_INTEGER, QueryPerformanceCounter() */
#endif

#ifdef _OPENMP
    #include <omp.h>
#endif


/* Abscissa of convergence and singularities of F(s) with their polar multiplicities: required */
unsigned int LTsings2 (unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0);

/* LT samples on Talbot's contour */
double complex *OMP_LTsamples_fun (unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS);

double *linspace(unsigned int Nv, double Vmin, double Vmax);


int main(int argc, char *argv[])
/** RUN AS:
            ./ex1.exe tol jFUN NTval NXval

        tol     tolerance (1e-6, 1e-8, 1e-10, 1e-12)
        jFUN    number of OMP function to be tested
                = 1 OMP_Talbot11 [coarse-grain parallelism for modified method]
                = 2 OMP_Talbot12 [ fine-grain  parallelism for modified method]
                = 3 OMP_Talbot13 [   nested    parallelism for modified method]
        NTval   number of t-values (5, 20, 120)
        NXval   number of x-values (5, 20, 120)
 **/
{
    /*   INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM        */
    double tol;                             /* tolerance      */
    unsigned int jFUN; /* selection of function to be tested  */
    unsigned int NTval, NXval;    /* number of t and x values */
    if (argc > 4)
    {
        tol   = (double)atof(argv[1]);       /* 1st argument to main */
        jFUN  = (unsigned int)atoi(argv[2]); /* 2nd argument to main */
        NTval = (unsigned int)atoi(argv[3]); /* 3rd argument to main */
        NXval = (unsigned int)atoi(argv[4]); /* 4th argument to main */
    }
    else
    {
        tol   = 1e-8;                        /* default values       */
        jFUN  = 1;
        NTval = 20;
        NXval = 20;
    }


    /* endpoints of intervals */
    double Xmin=0.0,    Xmax=5.0;

    double Tmin=100.0,   Tmax=500.0;


    char *STR = "OMP_LTsamples_fun";          /* OpenMP parallel LT samples */
    char *Talbot_fun[] = {"OMP_Talbot11_DE",  /* OpenMP coarse-grain parallelism for modified Talbot's method */
                          "OMP_Talbot12_DE",  /* OpenMP  fine-grain  parallelism for modified Talbot's method */
                          "OMP_Talbot13_DE"}; /* OpenMP    nested    parallelism for modified Talbot's method */

    unsigned int k; /* for-loop index */


    /* INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION */
    double sigma0;            /* abscissa of convergence                       */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    LTsings2(&Nsings, SINGS, MULT, &sigma0);

    char *LTS = "LT samples by means of a function";
    printf("%s", str);
    puts("\n\t   Ex. 2: output from ./2PAR/LTS1_fun/OMP_main_TIMES.c");
    printf("\n\t%s\n", LTS);
    printf("\n\t      t in [%g, %g],    x in [%g, %g],    tol=%e", Tmin,Tmax,Xmin,Xmax,tol);
    puts("\n====================================================================================\n");

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


    /* COMPUTE THE t AND x VALUES */
    double *t = linspace(NTval,Tmin,Tmax);
    double *x = linspace(NXval,Xmin,Xmax);


    /*   OUTPUT DATA                                       */
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot, *IFAIL; /* global and local error flags */

    /* ALLOCATE ALGORITHM'S OUTPUT ARRAYS */
    NUMft = (double *)malloc(NXval*NTval*sizeof(double));
    if (NUMft == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n"); exit(1);}

    IFAIL = (int *)malloc(NXval*NTval*sizeof(int));
    if (IFAIL == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n"); exit(1);}


    /* FOR TIMING [both Windows and Linux] */
    double TOT_TIME, MEAN_TIME;
    unsigned int NTIMES = 5, Nt;


    int THREADS, maxTHREADS;
    #ifdef _OPENMP
        printf("[OpenMP 4.0] omp_get_num_procs() = %d,   omp_get_max_threads() = %d,   omp_get_num_threads() = %d\n", omp_get_num_procs(),omp_get_max_threads(),omp_get_num_threads() );
        maxTHREADS = omp_get_num_procs(); /* OpenMP 4.0 */
    #else
        printf("OpenMP not defined!\n");
        maxTHREADS = 1;
    #endif

    char ch='%';
    printf("\n\t***   RESULTS OF OpenMP TALBOT SUITE DE function:   %s()   ***", Talbot_fun[jFUN-1]);
    printf("\n\t      %d t in [%g, %g],    %d x in [%g, %g],    input tolerance: tol = %5.2e", NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
    double CONLAM, CONSIG, CONNU; unsigned int NOPTS;
    COM_TalbotPAR (sigma0,(Tmin+Tmax)/2,tol,Nsings,SINGS,MULT, &CONLAM,&CONSIG,&CONNU,&NOPTS);
    NOPTS = COM_TalbotNcorr(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4);
    printf("\n\t      NOPTS = %u", NOPTS);
    printf("\n\t      user defined function to compute the LT samples:  %s()", STR);


    /* SELECTION OF FUNCTION TO BE TESTED */
    if ( jFUN == 3 )
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
        printf("\n\tMEAN ELAPSED TIME:\n");
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
        printf("\n\t      Modified Talbot's method [%s parallelism]", (jFUN == 1) ? "coarse-grain" : "fine-grain");
        printf("\n\tMEAN ELAPSED TIME:\n");
        switch (jFUN)
        {
            case 1 : /* (OpenMP coarse-grained parallel modified method of Talbot Suite DE)  */
                /* Loop on number of threads */
                for ( THREADS=1;  THREADS<=maxTHREADS;  THREADS++ )
                {
                    /* SET THE NUMBER OF OpenMP PARALLEL THREADS */
                    #ifdef _OPENMP
                        omp_set_dynamic(0); /* static scheduling */
                        omp_set_num_threads(THREADS);
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
                    /* PRINT RESULTS */
                    printf("\n\t\t%e\t%c number of threads = %2d", MEAN_TIME,ch,THREADS);
                } /* END for (THREADS) */
                puts("\n\n");
                break;

            case 2 : /* (OpenMP fine-grained parallel modified method of Talbot Suite DE)  */
                /* Loop on number of threads */
                for ( THREADS=1;  THREADS<=maxTHREADS;  THREADS++ )
                {
                    /* SET THE NUMBER OF OpenMP PARALLEL THREADS */
                    #ifdef _OPENMP
                        omp_set_dynamic(0); /* static scheduling */
                        omp_set_num_threads(THREADS);
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
                    /* PRINT RESULTS */
                    printf("\n\t\t%e\t%c number of threads = %2d", MEAN_TIME,ch,THREADS);
                } /* END for (THREADS) */
                puts("\n\n");
                break;

            default :
                fprintf(stderr, "\n*** No function corresponds to the value %d ***\n", jFUN);
                exit(1);
        } /* END switch jFUN */
    }

    return 0;
}


double *linspace(unsigned int Nv, double Vmin, double Vmax)
/** Compute a dynamic double array of equispace values
                    [sequential]
  */
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
