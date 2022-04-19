/********************        SEQ_main_ACCURACY.c       ********************
 *                                                                        *
 *                                                                        *
 *                        SEQUENTIAL TALBOT SUITE DE                      *
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
#include "../../../TalbotSuiteDE/FUN_DE/SEQ_Talbot_pack_DE.h" /* for Talbot Suite DE's sequential implementation */
#include "../../problem_string.h"


/*  Inverse LT function u(x,t): not required
    ========================================
    This function has been included to compare the numerical approximations
    to the corresponding true results
 */
double ILTfun2 (double x, double t); /* u(x,t) */

/* Abscissa of convergence and singularities of F(s) with their polar multiplicities: required */
unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0);

/* LT samples on Talbot's contour */
double complex *SEQ_LTsamples_ode (unsigned int NXval, double Xval[], unsigned int N, double complex S[], double tol);

double *rel_err (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol);
double *abs_err (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol);
double *linspace (unsigned int Nv, double Vmin, double Vmax);


int main(int argc, char *argv[])
/*  To call: ./ex5.exe 1e-6

    Parameter to the main function:
    tol: tolerance on error
*/
{
    char *LTS = "LT samples computed by solving ODE problems by means of ode.c";

    double tol;                             /* tolerance            */
    if (argc > 1)
        tol = (double)atof(argv[1]);        /* 1st argument to main */
    else
        tol = 1e-8;                         /* default value        */


    /*   INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM             */
    unsigned int NXval=9;       /* number of x values              */
    unsigned int NTval=5;       /* number of t values              */

    double Xmin=  0.0,  Xmax=  5.0; /* endpoints of the interval for x */
    double Tmin=100.0,  Tmax=500.0; /* endpoints of the interval for t */


    printf("%s", str);
    puts("\n\t   Ex. 2: output from ./1SEQ/LTS2_ode/SEQ_main_ACCURACY.c");
    printf("\n\t%s\n", LTS);
    printf("\n\t%d t in [%g, %g],    %d x in [%g, %g],    tol=%e", NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
    puts("\n====================================================================================\n");


    /* COMPUTE ERRORS */
    char ERtype = 'A'; /* 'A': absolute error
                          'R': relative error */
    double *ERR1, *ERR2;
    if (ERtype == 'A')
    {   ERR1 = abs_err (1,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol); /* SEQ_Talbot1_DE */
        ERR2 = abs_err (2,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol); /* SEQ_Talbot2_DE */
    }
    else /* relative error */
    {   ERR1 = rel_err (1,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol); /* SEQ_Talbot1_DE */
        ERR2 = rel_err (2,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol); /* SEQ_Talbot2_DE */
    }


    /* DISPLAY ERROR MATRICES */
    unsigned int h, k; /* for-loop index */
    char ch='%';
    (ERtype == 'A') ? printf("\nABS") :  printf("\nREL") ;
    printf("ERR1 = [ %c   Tval(1)        Tval(2)        Tval(3) ...\n", ch);
    for (h=0; h<NXval; h++)
    {   printf("          ");
        for (k=0; k<NTval; k++)    printf("  %e", ERR1[h*NTval+k]);
        printf("    %c Xval(%2d)\n", ch,h+1);
    }
    printf(  "          ];\n");

    (ERtype == 'A') ? printf("\nABS") :  printf("\nREL") ;
    printf("ERR2 = [ %c   Tval(1)        Tval(2)        Tval(3) ...\n", ch);
    for (h=0; h<NXval; h++)
    {   printf("          ");
        for (k=0; k<NTval; k++)    printf("  %e", ERR2[h*NTval+k]);
        printf("    %c Xval(%2d)\n", ch,h+1);
    }
    printf(  "          ];\n");

    return 0;
}


double *rel_err (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol)
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
        case 1 : /* (modified Talbot's method */
            IFAIL_tot = SEQ_Talbot1_DE ( SEQ_LTsamples_ode, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax);
            break;

        case 2 : /* classical Talbot's method */
            IFAIL_tot = SEQ_Talbot2_DE ( SEQ_LTsamples_ode, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT);
            break;

        default :
            fprintf(stderr, "\n*** No function corresponds to the value %d ***", jFUN);
            exit(1);

    } /* END switch jFUN */

    if ( IFAIL_tot != 0 )
        printf("Total error indicator: IFAIL_tot = %d", IFAIL_tot);

    /* RELATIVE ERRORS */
    double ft;
    double *RELERR = (double *)malloc(NXval*NTval*sizeof(double));
    if (RELERR == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF RELERR IS FAILED. ***\n"); exit(1);}

    for (h=0; h<NXval; h++)
        for (k=0; k<NTval; k++)
        {   ft = ILTfun2( x[h], t[k] );
            RELERR[h*NTval+k] = fabs( ft - NUMft[h*NTval + k] )/fabs(ft);
        }

    free(t); free(x); free(NUMft); free(IFAIL);
    return RELERR;
}


double *abs_err (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol)
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
        case 1 : /* (modified Talbot's method */
            IFAIL_tot = SEQ_Talbot1_DE ( SEQ_LTsamples_ode, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax);
            break;

        case 2 : /* classical Talbot's method */
            IFAIL_tot = SEQ_Talbot2_DE ( SEQ_LTsamples_ode, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT);
            break;

        default :
            fprintf(stderr, "\n*** No function corresponds to the value %d ***", jFUN);
            exit(1);

    } /* END switch jFUN */

    if ( IFAIL_tot != 0 )
        printf("Total error indicator: IFAIL_tot = %d", IFAIL_tot);

    /* ABSOLUTE ERRORS */
    double ft;
    double *ABSERR = (double *)malloc(NXval*NTval*sizeof(double));
    if (ABSERR == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF ABSERR IS FAILED. ***\n"); exit(1);}

    for (h=0; h<NXval; h++)
        for (k=0; k<NTval; k++)
        {   ft = ILTfun2( x[h], t[k] );
            ABSERR[h*NTval+k] = fabs( ft - NUMft[h*NTval + k] );
        }

    free(t); free(x); free(NUMft); free(IFAIL);
    return ABSERR;
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
