/********************        SEQ_main_ACCURACY.c       ********************
 *                                                                        *
 *                                                                        *
 *                        SEQUENTIAL TALBOT SUITE DE                      *
 *                                                                        *
 *                      DRIVER PROGRAM FOR EXAMPLE 4a                     *
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
#include <math.h>
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

double *abs_rel_err (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol);
double *linspace (unsigned int Nv, double Vmin, double Vmax);


int main(int argc, char *argv[])
/*  ./ex4_acc.exe 1e-6

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
/*    unsigned int NXval=10; */
/*    unsigned int NTval=10; */

    /* ENDPOINTS */
    double  Tmin=100.0,   Tmax=500.0;
    double  Xmin=0.0,     Xmax=1.0;


    printf("%s", str);
    puts("\n\t   Ex. 4a: output from ./1SEQ/LTS2_ode/SEQ_main_ACCURACY.c");
    printf("\n\t%s\n", LTS);
    printf("\n\t%d t in [%g, %g],    %d x in [%g, %g],    tol=%e", NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
    puts("\n====================================================================================\n");


    /* COMPUTE ABSOLUTE AND RELATIVE ERRORS */
    double *ERROR1, *ERROR2;
    ERROR1 = abs_rel_err (1,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol); /* SEQ_Talbot1_DE */
    ERROR2 = abs_rel_err (2,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol); /* SEQ_Talbot2_DE */


    /* DISPLAY ERROR MATRICES */
    char ch='%';
    unsigned int h, k; /* for-loop index */

    /* ABSOLUTE ERRORS */
/*
    printf("\nABSERR1 = [ %c  Tval(1)    Tval(2)    ...    Tval(%d)\n", ch,NTval);
    for (h=0; h<NXval; h++)
    {   printf("          ");
        for (k=0; k<NTval; k++)    printf("  %e", ERROR1[h*NTval+k]);
        printf("    %c Xval(%2d)\n", ch,h+1);
    }
    printf(  "          ];\n");
    printf("\nABSERR2 = [ %c  Tval(1)    Tval(2)    ...    Tval(%d)\n", ch,NTval);
    for (h=0; h<NXval; h++)
    {   printf("          ");
        for (k=0; k<NTval; k++)    printf("  %e", ERROR2[h*NTval+k]);
        printf("    %c Xval(%2d)\n", ch,h+1);
    }
    printf(  "          ];\n");
*/

    /* RELATIVE ERRORS */
    printf("\nRELERR1 = [ %c  Tval(1)    Tval(2)    ...    Tval(%d)\n", ch,NTval);
    for (h=0; h<NXval; h++)
    {   printf("          ");
        for (k=0; k<NTval; k++)    printf("  %e", ERROR1[(NXval+h)*NTval+k]);
        printf("    %c Xval(%2d)\n", ch,h+1);
    }
    printf(  "          ];\n");
    printf("\nRELERR2 = [ %c  Tval(1)    Tval(2)    ...    Tval(%d)\n", ch,NTval);
    for (h=0; h<NXval; h++)
    {   printf("          ");
        for (k=0; k<NTval; k++)    printf("  %e", ERROR2[(NXval+h)*NTval+k]);
        printf("    %c Xval(%2d)\n", ch,h+1);
    }
    printf(  "          ];\n");

    printf("\n\t      %d t in [%g, %g],    %d x in [%g, %g],    tol = %e\n", NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
    return 0;
}


double *abs_rel_err (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol)
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
    double sigma0;             /* abscissa of convergence                       */
    unsigned int Nsings;       /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10];  /* array for singularities                       */
    unsigned int MULT[10];     /* array for polar multiplicities                */
    LTsings2(&Nsings, SINGS, MULT, &sigma0);


    /* TALBOT SUITE DE user level functions  */
    switch (jFUN)
    {
        case 1 : /* (modified Talbot's method */
            IFAIL_tot = SEQ_Talbot1_DE ( SEQ_LTsamples_ode, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,
                                         Nsings,SINGS,MULT,Tmin,Tmax );
            break;

        case 2 : /* classical Talbot's method */
            IFAIL_tot = SEQ_Talbot2_DE ( SEQ_LTsamples_ode, sigma0,NXval,x,NTval,t,tol,NUMft,IFAIL,
                                         Nsings,SINGS,MULT );
            break;

        default :
            fprintf(stderr, "\n*** No function corresponds to the value %d ***", jFUN);
            exit(1);

    } /* END switch jFUN */

    if ( IFAIL_tot != 0 )
        printf("Total error indicator: IFAIL_tot = %d", IFAIL_tot);

    /*  ABSOLUTE AND RELATIVE ERRORS
        ERROR is a matrix of size (2*NXval,NTval) such that:
                        |                   |
                        |  absolute errors  |
                        |                   |
                ERROR = ---------------------
                        |                   |
                        |  relative errors  |
                        |                   |
        * the upper half matrix contains the absolute error
                ERROR[h*NTval+k]: absolute error in u(x[h],t[k])
        * the lower half matrix contains the relative error
                ERROR[(NXval+h)*NTval+k]: relative error in u(x[h],t[k])
     */
    double ft;
    double *ERROR = (double *)malloc(2*NXval*NTval*sizeof(double));
    if (ERROR == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF ERROR matrix IS FAILED. ***\n"); exit(1);}

    for (h=0; h<NXval; h++)
        for (k=0; k<NTval; k++)
        {   ft = ILTfun2( x[h], t[k] );
            ERROR[h*NTval+k] = fabs( ft - NUMft[h*NTval + k] );   /* absolute error */
            ERROR[(NXval+h)*NTval+k] = ERROR[h*NTval+k]/fabs(ft); /* relative error */
        }

    free(t); free(x); free(NUMft); free(IFAIL);
    return ERROR;
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
