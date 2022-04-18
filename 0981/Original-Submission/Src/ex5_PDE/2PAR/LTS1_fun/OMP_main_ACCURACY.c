/********************        OMP_main_ACCURACY.c       ********************
 *                                                                        *
 *                                                                        *
 *                  OpenMP-BASED PARALLEL TALBOT SUITE DE                 *
 *                                                                        *
 *                       DRIVER PROGRAM FOR EXAMPLE 5                     *
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
#include "../../../TalbotSuiteDE/FUN_DE/OMP_Talbot_pack_DE.h" /* for Talbot Suite DE's OpenMP implementation */
#include "../../problem_string.h"

#ifdef _OPENMP
    #include <omp.h>
#endif


/*  Inverse LT function u(x,t): not required
    ========================================
    This function has been included to compare the numerical approximations
    to the corresponding true results
 */
double ILTfun2 (double x, double y, double t); /* u(x,y,t) */

/* Abscissa of convergence and singularities of F(s) with their polar multiplicities: required */
unsigned int LTsings2(unsigned int *Nsings, double complex *SINGS, unsigned int *MULT, double *sigma0);

/* LT samples on Talbot's contour */
double complex *OMP_LTsamples_fun (unsigned int NXYval, double XY[], unsigned int N, double complex S[], double tol, int THREADS);

double *internalMeshPts(unsigned int NXYval, double XYmin, double XYmax);

double *rel_err (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXval, double Xmin, double Xmax, double tol, int thrds1, int thrds2);
void display_err(char *STR, char FUNtype, unsigned int NXval, unsigned int NTval, double *ERROR, int thrds1, int thrds2);
double *linspace (unsigned int Nv, double Vmin, double Vmax);


int main(int argc, char *argv[])
/** RUN AS:
            ./ex5_acc.exe tol thrds1
    OR AS:
            ./ex5_acc.exe tol thrds1 thrds2     [for nested parallelism]
**/
{
    double tol;                      /* tolerance                */
    int thrds1, thrds2;              /* number of OpenMP threads */
    switch ( argc ) {

    case 4 :
        /* ./ex5_acc.exe tol threads1 threads2 [nested parallelism] */
        tol = (double)atof(argv[1]); /* 1st argument to main */
        thrds1 = (int)atoi(argv[2]); /* 2nd argument to main */
        thrds2 = (int)atoi(argv[3]); /* 3rd argument to main */
        break;

    case 3 :
        /* ./ex5_acc.exe tol threads */
        tol = (double)atof(argv[1]); /* 1st argument to main */
        thrds1 = (int)atoi(argv[2]); /* 2nd argument to main */
        thrds2 = 1;
        break;

    default :
        tol = 1e-8;
        thrds1 = 1;
        thrds2 = 1;
    }


    /*  SET THE NUMBER OF OpenMP PARALLEL THREADS */
    #ifdef _OPENMP
        omp_set_dynamic(0); /* static scheduling */
        omp_set_num_threads(thrds1*thrds2);
    #endif

    /*   INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM                 */
    unsigned int NXYval=9;         /* number of x,y values            */
    unsigned int NTval =5;         /* number of t values              */


    /* endpoints */
    double XYmin=0.0,     XYmax=1.0;
    double Tmin=100.0,    Tmax=500.0;


    char *LTS = "LT samples by a function";
    printf("%s", str);
    puts("\n\t   Ex. 5: output from ./2PAR/LTS1_fun/OMP_main_ACCURACY.c");
    printf("\n\t%s\n", LTS);
    printf("\n\t      %d t in [%g, %g],    %d x,y in [%g, %g],    tol=%e", NTval,Tmin,Tmax,NXYval,XYmin,XYmax,tol);
    printf("\n====================================================================================\n\n");


    /*  COMPUTE ERRORS in u(x,y,t)
            u(j,k) = u(X(j),Y(j),t(k))
     */
    unsigned int Nrows = (NXYval-2)*(NXYval-2);
    double *RELERR1, *RELERR2, *RELERR3;
    RELERR1 = rel_err (1,NTval,Tmin,Tmax,NXYval,XYmin,XYmax,tol,thrds1,thrds2);
    RELERR2 = rel_err (2,NTval,Tmin,Tmax,NXYval,XYmin,XYmax,tol,thrds1,thrds2);
    RELERR3 = rel_err (3,NTval,Tmin,Tmax,NXYval,XYmin,XYmax,tol,thrds1,thrds2);

    /* DISPLAY ERROR MATRICES */
    display_err("REL",'1',Nrows,NTval,RELERR1,thrds1,thrds2);
    display_err("REL",'2',Nrows,NTval,RELERR2,thrds1,thrds2);
    display_err("REL",'3',Nrows,NTval,RELERR3,thrds1,thrds2);

    free(RELERR1); free(RELERR2); free(RELERR3);

    printf("\n\t      %d t in [%g, %g],    %d x,y in [%g, %g],    tol=%e\n\n", NTval,Tmin,Tmax,NXYval,XYmin,XYmax,tol);
    return 0;
}


double *rel_err (unsigned int jFUN, unsigned int NTval, double Tmin, double Tmax, unsigned int NXYval, double XYmin, double XYmax, double tol, int thrds1, int thrds2)
{
    /* CARTESIAN COORDINATES OF INTERNAL MESH POINTS */
    double *XY = internalMeshPts(NXYval,XYmin,XYmax);

    unsigned int h, k; /* for loop index */
    unsigned int Nrows = (NXYval-2)*(NXYval-2); /* number of internal mesh points */

    /*  OUTPUT DATA
        NUMft, IFAIL: output matrices of size of ((NXYval-2)^2, NTval))
    */
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot, *IFAIL; /* global and local error flags */
    NUMft = (double *)malloc(Nrows*NTval*sizeof(double));
    if (NUMft == NULL)
    {   fprintf(stderr,"\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n");
        exit(1);
    }

    IFAIL = (int *)malloc(Nrows*NTval*sizeof(int));
    if (IFAIL == NULL)
    {   fprintf(stderr,"\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n");
        exit(1);
    }


    /* INPUT DATA ABOUT LAPLACE TRANSFORM */
    double sigma0;            /* abscissa of convergence                       */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    LTsings2(&Nsings, SINGS, MULT, &sigma0);


    /* INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM */
    double *t = linspace(NTval,Tmin,Tmax);


    /* TALBOT SUITE DE user level functions  */
    switch (jFUN)
    {
        case 1 : /* modified Talbot's method coarse grain parallelism */
            IFAIL_tot = OMP_Talbot11_DE ( OMP_LTsamples_fun, sigma0,Nrows,XY,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax,thrds1);
            break;

        case 2 : /* modified Talbot's method fine-grain parallelism */
            IFAIL_tot = OMP_Talbot12_DE ( OMP_LTsamples_fun, sigma0,Nrows,XY,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax,thrds1);
            break;

        case 3 : /* modified Talbot's method hybrid fine/coarse-grain parallelism */
            #ifdef _OPENMP
                printf("\nCHECK: NESTED PARALLELISM is %s\n", omp_get_nested() ? "supported" : "not supported");
                if ( !omp_get_nested() )
                {   omp_set_nested(1); /* true */
                    printf("                          Now it is %s\n", omp_get_nested() ? "supported" : "not supported");
                }
            #endif
            IFAIL_tot = OMP_Talbot13_DE ( OMP_LTsamples_fun, sigma0,Nrows,XY,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax,thrds1,thrds2);
            break;

        default :
        {   fprintf(stderr, "\n*** No function corresponds to the value %d ***", jFUN);
            exit(1);
        }

    } /* END switch jFUN */

    if ( IFAIL_tot != 0 )
        printf("Total error indicator: IFAIL_tot = %d\n", IFAIL_tot);


    /* RELATIVE ERRORS */
    double ft;
    double *RELERR = (double *)malloc(Nrows*NTval*sizeof(double));
    if (RELERR == NULL)
    {   fprintf(stderr,"\n***   ERROR: DYNAMIC ALLOCATION OF RELERR IS FAILED. ***\n");
        exit(1);
    }

    for (h=0; h<Nrows; h++)
        for (k=0; k<NTval; k++)
        {   ft = ILTfun2( XY[h], XY[Nrows+h], t[k] ); /* col-wise XY(h,j) = XY[j*Nrows+h]) */
            RELERR[h*NTval + k] = fabs( ft - NUMft[h*NTval + k] )/fabs(ft);
        }

    free(t); free(XY); free(NUMft); free(IFAIL);
    return RELERR;
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
    {   printf("%2d        ", h+1);
        for (k=0; k<NTval; k++)    printf("  %e", ERROR[h*NTval+k]);
        printf("\n");
    }
    printf(  "          ];\n");
}


double *linspace(unsigned int Nv, double Vmin, double Vmax)
/** Compute a dynamic double array of equispace values **/
{
    double *v = (double *)malloc(Nv*sizeof(double));
    if (v == NULL)
    {   fprintf(stderr,"\n***   ERROR IN linspace(): DYNAMIC ALLOCATION OF THE ARRAY IS FAILED. ***\n");
        exit(1);
    }

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
