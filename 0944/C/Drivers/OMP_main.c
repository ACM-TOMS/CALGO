/************************        OMP_main.c       *************************
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *             DRIVER PROGRAM FOR THE OMP-BASED IMPLEMENTATION            *
 *                             OF TALBOT SUITE                            *
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *  >>>>>>>>>>>>>       VERSION 1.0   Sept 13th, 2012         <<<<<<<<<<  *
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *           AUTHORS: Laura Antonelli (1), Stefania Corsaro (2-1),        *
 *                    Zelda Marino (2), Mariarosaria Rizzardi (3)         *
 *                                                                        *
 *                   (1) ICAR  - National Research Council of Italy       *
 *                                                                        *
 *                   (2) DSMRE - "Parthenope" University, Naples (Italy)  *
 *                                                                        *
 *                   (3) DSA   - "Parthenope" University, Naples (Italy)  *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 * REFERENCES                                                             *
 * ==========                                                             *
 * Antonelli L., Corsaro S.,                                              *
 * Marino Z., Rizzardi M. - "Talbot Suite: a parallel software collection *
 *                    for the numerical inversion of Laplace Transforms". *
 *                    ACM Trans. Math. Softw., vol.##, no.#, month year,  *
 *                    pp. ##-##.                                          *
 *                                                                        *
 **************************************************************************/

/* Required GNU gcc compiler option:    -std=gnu99 [for clock_gettime()]
   Required         linker   option:    -lrt       [for clock_gettime()]
*/

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "../SRC/main_include.h"
#include "../SRC/OMP_Talbot_pack.h" /* for Talbot Suite's OpenMP implementation */
#include "../SRC/COM_Talbot_pack.h" /* for Talbot Suite's shared functions      */
#include <time.h>                   /* for timing                               */


/*  Laplace Transform function (LTfun)
    [F24(s) from test functions]:
        F(s) = F(s) = s/(s^2+9)^2

    2 polar singularities (double poles):
        s_j = {+/-3*i}

    abscissa of convergence:
        sigma0 = 0

    Inverse Laplace Transform function (ILTfun):
        f(t) = t*sin(3*t)/6
*/
double complex LTfun(double complex s)
/* Laplace Transform function F(s) */
{   double complex LapVal = s*s+9;
    return s/(LapVal*LapVal);
}

double ILTfun(double t)
/* Inverse Laplace Transform function f(t) */
{   return t*sin(3.0*t)/6;
}

unsigned int LTsings(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0)
/* Abscissa of convergence and singularities of F(s) with their polar multiplicities */
{   unsigned int NsingsTOT = 2; /* total number of singularities s_j */
    *Nsings = 1;                /* number of singularities s_j with Im(s_j) ge 0 */
    *sigma0 = 0;                /* abscissa of convergence */

    SINGS[0]=0.0+3.0*I; SINGS[1]=0.0-3.0*I;
    MULT[0]=2; MULT[1]=2;

    return NsingsTOT;
}



int main(int argc, char *argv[])
/* Parameters to the main function:
    1) jFUN   : switch between OMP_Talbot1 (1) and OMP_Talbot2 (2)
    2) THREADS: number of threads in parallel sections
*/
{
    int k; /* for-loop index */


    /* SELECTION OF TALBOT SUITE FUNCTION */
    char *Talbot_fun[] = {"OMP_Talbot1",  /* coarse grain parallelism  */
                          "OMP_Talbot2"}; /*  fine  grain parallelism  */
    unsigned int jFUN; /* switch between OMP_Talbot1 (1)
                                     and OMP_Talbot2 (2) */
    if (argc >= 2)
        jFUN = (unsigned int)atoi(argv[1]); /* 1st argument to main */
    else
        jFUN = 1;                           /* default value        */


    /* SET OpenMP ENVIRONMENT */
    int THREADS; /* number of threads in parallel sections */
    #ifdef _OPENMP
        omp_set_dynamic(0); /* static scheduling */
        /* SET THE NUMBER OF THREADS */
        if (argc == 3)
            THREADS = (int)atoi(argv[2]); /* 2nd argument to main */
        else
            THREADS = 1;                  /* default value        */
        omp_set_num_threads(THREADS);
    #else
        THREADS = 1;                      /* run as sequential    */
    #endif


    /* INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION */
    char *LTstring  = "F(s) = s/(s^2+9)^2";
    char *ILTstring = "f(t) = t*sin(3*t)/6";
    double sigma0;            /* abscissa of convergence                       */
    unsigned int NsingsTOT;   /* total number of singularities s_j             */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    NsingsTOT = LTsings(&Nsings, SINGS, MULT, &sigma0);

    printf("\nLaplace Transform function:\t\t%s\tabscissa of convergence: sigma0 = %f\n", LTstring,sigma0);
    printf("\nInv. Laplace Transform function:\t%s\n\n", ILTstring);
    (NsingsTOT == 1  ?  printf("   singularity    and    multiplicity:\n")
	                 : printf("   singularities    and    multiplicities:\n"));
    for (k=0; k<NsingsTOT; k++)
        (fabs(cimag(SINGS[k]))>0 ? printf("\ns(%2d) =%+6.2f%+6.2f * I\t\tmult = %d", k+1,creal(SINGS[k]),cimag(SINGS[k]),MULT[k])
		                         : printf("\ns(%2d) =%+6.2f\t\t\tmult = %d", k+1,creal(SINGS[k]),MULT[k]));
    puts("\n");


    /* INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM                 */
    unsigned int NTval=20;     /* number of t values                 */
    double Tmin=100, Tmax=500; /* endpoints of the interval for t    */
    double tol = 1e-12;        /* tolerance                          */
    double *t;                 /* t array                            */

    /* OUTPUT DATA */
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot, *IFAIL; /* global and local error flags */


    /* COMPUTE THE t VALUES */
    t = (double *)malloc(NTval*sizeof(double));
    if (t == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF t IS FAILED. ***\n"); exit(1);}
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
    struct timespec ts1, ts2; double T0, T1, TOT_TIME; /* for timing */


    /* CALL THE TALBOT SUITE FUNCTION  */
    switch (jFUN)
    {
        case 1 : /* coarse grain parallelism (modified Talbot's method)  */
            clock_gettime(CLOCK_REALTIME, &ts1);
                IFAIL_tot = OMP_Talbot1 (LTfun,sigma0,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax);
            clock_gettime(CLOCK_REALTIME, &ts2);
            T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec      */
            T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
            TOT_TIME = T0 + T1*1.0e-9;
            break;

        case 2 : /*  fine  grain parallelism (classical Talbot's method) */
            clock_gettime(CLOCK_REALTIME, &ts1);
                IFAIL_tot = OMP_Talbot2 (LTfun,sigma0,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT);
            clock_gettime(CLOCK_REALTIME, &ts2);
            T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec      */
            T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
            TOT_TIME = T0 + T1*1.0e-9;
            break;

        default :
            fprintf(stderr, "\n*** No function corresponds to the value %d ***", jFUN);
            exit(1);

    } /* END switch jFUN */


    /* PRINT RESULTS */
    printf("\nInverting the LT fun at NTval=%d values of t in [%5.2f, %5.2f]\n", NTval,Tmin,Tmax);
    printf("\n\t***   RESULTS OF OMP PARALLEL TALBOT SUITE [function = %s()]   ***\n", Talbot_fun[jFUN-1]);
    printf("\n\t\tinput tolerance: tol = %5.2e\t\tthreads number: %d\telapsed time = %e\n", tol,THREADS,TOT_TIME);
    if (NTval < 50)
    {
        if (!IFAIL_tot)
        {
            printf("\n\n       T       F EXACT          F APPROX        ABS ERR         REL ERR      TYPE    IFAIL_tot = %d (no local error)\n", IFAIL_tot);
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
            puts("\n\n       T       F EXACT          F APPROX        ABS ERR         REL ERR      TYPE    IFAIL");
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
    return 0;
}
