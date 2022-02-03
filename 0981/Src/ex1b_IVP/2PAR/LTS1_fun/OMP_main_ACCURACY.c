/********************        OMP_main_ACCURACY.c       ********************
 *                                                                        *
 *                                                                        *
 *                  OpenMP-BASED PARALLEL TALBOT SUITE DE                 *
 *                                                                        *
 *                      DRIVER PROGRAM FOR EXAMPLE 1b                     *
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
#include "../../../TalbotSuiteDE/FUN_DE/OMP_Talbot_pack_DE.h" /* for Talbot Suite DE's OpenMP parallel implementation */
#include "../../problem_string.h"

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
double complex *OMP_LTsamples_fun (unsigned int NXval, double Xval[], unsigned int N, double complex S[], double tol, int THREADS);

void errors3 (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax,
              double tol, double *ABSERR, double *RELERR, double *PS_ERR, int thrds1, int thrds2);
void display_err(char *STR, char FUNtype, unsigned int NXval, unsigned int NTval, double *ERROR, int thrds1, int thrds2);
double *linspace (unsigned int Nv, double Vmin, double Vmax);


int main(int argc, char *argv[])
/** RUN AS:
            ./ex1_acc.exe tol thrds1
    OR AS:
            ./ex1_acc.exe tol thrds1 thrds2     [for nested parallelism]
**/
{
    double tol;                      /* tolerance                */
    int thrds1, thrds2;              /* number of OpenMP threads */
    switch ( argc ) {

    case 4 :
        /* ./ex1_acc.exe tol threads1 threads2 [nested parallelism] */
        tol = (double)atof(argv[1]); /* 1st argument to main */
        thrds1 = (int)atoi(argv[2]); /* 2nd argument to main */
        thrds2 = (int)atoi(argv[3]); /* 3rd argument to main */
        break;

    case 3 :
        /* ./ex1_acc.exe tol threads */
        tol = (double)atof(argv[1]); /* 1st argument to main */
        thrds1 = (int)atoi(argv[2]); /* 2nd argument to main */
        thrds2 = 1;
        break;

    default :
        tol = 1e-8;
        thrds1 = 1;
        thrds2 = 1;
    }


    /* SET THE NUMBER OF OpenMP PARALLEL THREADS */
    #ifdef _OPENMP
        omp_set_dynamic(0); /* static scheduling */
        omp_set_num_threads(thrds1*thrds2);
    #endif


    /*   INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM                 */
    unsigned int NXval=9;           /* number of x values              */
    double Xmin= 10.0,  Xmax= 20.0; /* endpoints of the interval for x */
    unsigned int NTval=5;           /* number of t values              */
    double Tmin=100, Tmax=500;      /* endpoints of the interval for t */


    char *LTS = "LT samples by means of a function";
    printf("%s", str);
    puts("\n\t   Ex. 1b: output from ./2PAR/LTS1_fun/OMP_main_ACCURACY.c");
    printf("\n\t%s\n", LTS);
    printf("\n\t      %d t in [%g, %g],    %d x in [%g, %g],    tol=%e", NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
    puts("\n====================================================================================\n");


    /* COMPUTE ERRORS */
    double *ABSERR1 = (double *)malloc(NXval*NTval*sizeof(double)),
           *ABSERR2 = (double *)malloc(NXval*NTval*sizeof(double)),
           *ABSERR3 = (double *)malloc(NXval*NTval*sizeof(double));
    if ( ABSERR1 == NULL  ||  ABSERR2 == NULL  ||  ABSERR3 == NULL )
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF ABSERR IS FAILED. ***\n"); exit(1);}

    double *RELERR1 = (double *)malloc(NXval*NTval*sizeof(double)),
           *RELERR2 = (double *)malloc(NXval*NTval*sizeof(double)),
           *RELERR3 = (double *)malloc(NXval*NTval*sizeof(double));
    if ( RELERR1 == NULL  ||  RELERR2 == NULL  ||  RELERR3 == NULL )
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF RELERR IS FAILED. ***\n"); exit(1);}

    double *PS_ERR1 = (double *)malloc(NXval*NTval*sizeof(double)),
           *PS_ERR2 = (double *)malloc(NXval*NTval*sizeof(double)),
           *PS_ERR3 = (double *)malloc(NXval*NTval*sizeof(double));
    if ( PS_ERR1 == NULL  ||  PS_ERR2 == NULL  ||  PS_ERR3 == NULL )
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF PS_ERR IS FAILED. ***\n"); exit(1);}

    errors3 (1,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol,ABSERR1,RELERR1,PS_ERR1,thrds1,thrds2);
    errors3 (2,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol,ABSERR2,RELERR2,PS_ERR2,thrds1,thrds2);
    errors3 (3,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol,ABSERR3,RELERR3,PS_ERR3,thrds1,thrds2);

    /* DISPLAY ERROR MATRICES */
    display_err("REL",'1',NXval,NTval,RELERR1,thrds1,thrds2);
    display_err("REL",'2',NXval,NTval,RELERR2,thrds1,thrds2);
    display_err("REL",'3',NXval,NTval,RELERR3,thrds1,thrds2);
    display_err("ABS",'1',NXval,NTval,ABSERR1,thrds1,thrds2);
    display_err("ABS",'2',NXval,NTval,ABSERR2,thrds1,thrds2);
    display_err("ABS",'3',NXval,NTval,ABSERR3,thrds1,thrds2);
    display_err("PS_",'1',NXval,NTval,PS_ERR1,thrds1,thrds2);
    display_err("PS_",'2',NXval,NTval,PS_ERR2,thrds1,thrds2);
    display_err("PS_",'3',NXval,NTval,PS_ERR3,thrds1,thrds2);

    free(ABSERR1); free(ABSERR2); free(ABSERR3);
    free(RELERR1); free(RELERR2); free(RELERR3);
    free(PS_ERR1); free(PS_ERR2); free(PS_ERR3);

    printf("\n\t      %d t in [%g, %g],    %d x in [%g, %g],    tol = %e\n", NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
    return 0;
}


void display_err(char *STR, char FUNtype, unsigned int NXval, unsigned int NTval, double *ERROR, int thrds1, int thrds2)
/** STR = "ABS" or "REL" or "PS_"
    FUNtype = '1', '2' or '3'
    ERROR: row-wise matrix of size (NXval,NTval)
 */
{
    char ch = '%';
    unsigned int h, k; /* for-loop index */

    printf("\n%sERR%c = [", STR,FUNtype);
    #ifdef _OPENMP
        switch ( FUNtype )
        {
        case '1' :
            printf(" %c COARSE-GRAIN OMP PARALLELISM with %d threads for modified Talbot's method\n", ch,thrds1 );
            break;
        case '2' :
            printf(" %c FINE-GRAIN OMP PARALLELISM with %d threads for modified Talbot's method\n", ch,thrds1 );
            break;
        case '3' :
            printf(" %c NESTED COARSE/FINE-GRAIN OMP PARALLELISM with %d (outer), %d (inner) threads for modified Talbot's method\n", ch,thrds1,thrds2 );
        }
    #else
        printf("\n");
    #endif
    for (h=0; h<NXval; h++)
    {   printf("          ");
        for (k=0; k<NTval; k++)    printf("  %e", ERROR[h*NTval+k]);
        printf("\n");
    }
    printf(  "          ];\n");
}


void errors3 (unsigned int jFUN,unsigned int NTval,double Tmin,double Tmax,unsigned int NXval,double Xmin,double Xmax,double tol,
              double *ABSERR,double *RELERR,double *PS_ERR,int thrds1,int thrds2)
{
    /* INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM */
    double *t = linspace(NTval,Tmin,Tmax);

    /* INPUT DATA ABOUT ODE, COMING FROM APPLYING THE LAPLACE SUBSTITUTION METHOD TO THE PDE */
    double *x = linspace(NXval,Xmin,Xmax);


    /* OUTPUT DATA */
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot, *IFAIL; /* global and local error flags */
    NUMft = (double *)malloc(NXval*NTval*sizeof(double));
    if (NUMft == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n"); exit(1);}

    IFAIL = (int *)malloc(NXval*NTval*sizeof(int));
    if (IFAIL == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n"); exit(1);}


    /* INPUT DATA ABOUT LAPLACE TRANSFORM */
    double sigma0;            /* abscissa of convergence                       */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    LTsings2(&Nsings, SINGS, MULT, &sigma0);


    /* TALBOT SUITE DE user level functions  */
    switch (jFUN)
    {
        case 1 : /* (modified Talbot's method coarse-grain parallelism */
            IFAIL_tot = OMP_Talbot11_DE ( OMP_LTsamples_fun, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax,thrds1);
            break;

        case 2 : /* modified Talbot's method fine-grain parallelism */
            IFAIL_tot = OMP_Talbot12_DE ( OMP_LTsamples_fun, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax,thrds1);
            break;

        case 3 : /* modified Talbot's method hybrid fine/coarse-grain parallelism */
            #ifdef _OPENMP
                printf("\nCHECK: NESTED PARALLELISM is %s\n", omp_get_nested() ? "supported" : "not supported");
                if ( !omp_get_nested() )
                {   omp_set_nested(1); /* true */
                    printf("                          Now it is %s\n", omp_get_nested() ? "supported" : "not supported");
                }
            #endif
            IFAIL_tot = OMP_Talbot13_DE ( OMP_LTsamples_fun, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax,thrds1,thrds2);
            break;

        default :
            fprintf(stderr, "\n*** No function corresponds to the value %d ***\n", jFUN);
            exit(1);

    } /* END switch jFUN */

    if ( IFAIL_tot != 0 )
        printf("Total error indicator: IFAIL_tot = %d", IFAIL_tot);


    /* COMPUTE ABSOLUTE, RELATIVE AND PSEUDO ERRORS */
    double ft;
    unsigned int h, k;
    for (h=0; h<NXval; h++)
        for (k=0; k<NTval; k++)
        {   ft = ILTfun2( x[h], t[k] ); /* true value */
            ABSERR[h*NTval+k] = fabs( ft - NUMft[h*NTval + k] );
            RELERR[h*NTval+k] = ABSERR[h*NTval+k]/fabs(ft);
            if (fabs(ft) > 1)
                PS_ERR[h*NTval+k] = RELERR[h*NTval+k];
            else
                PS_ERR[h*NTval+k] = ABSERR[h*NTval+k];
        }

    free(t); free(x); free(NUMft); free(IFAIL);
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
