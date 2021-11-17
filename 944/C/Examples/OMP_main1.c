/************************       OMP_main1.c       *************************
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *                 SAMPLE MAIN for OMP_Talbot1() function                 *
 *            IN THE OMP-BASED IMPLEMENTATION OF TALBOT SUITE             *
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

/* Required GNU gcc compiler option:    -std=c99   [for complex type]
*/

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "../SRC/main_include.h"
#include "../SRC/OMP_Talbot_pack.h" /* for Talbot Suite's OpenMP implementation */
#include "../SRC/COM_Talbot_pack.h" /* for Talbot Suite's shared functions      */

/*  Laplace Transform function:
        F(s) = s^3/(s^4+4)

    4 polar singularities:
        s_j = {-1+i,+1+i,-1+i,+1-i}

    abscissa of convergence:
        sigma0 = +1

    Inverse Laplace Transform function:
        f(t) = cos(t)*cosh(t)
*/
double complex LTfun(double complex s)
/* Laplace Transform function F(s) = s^3/(s^4 + 4) */
{   double complex LapVal = s*s*s;
    return LapVal/(LapVal*s+4.0);
}

int main(int argc, char *argv[])
/* Parameter to the main function:
    1) THREADS: number of threads in parallel sections
*/
{
    int k; /* for-loop index */


    /* SET OpenMP ENVIRONMENT */
    int THREADS; /* number of threads in parallel sections */
    #ifdef _OPENMP
        omp_set_dynamic(0); /* static scheduling */
        /* SET THE NUMBER OF THREADS */
        if (argc == 2)
            THREADS = (int)atoi(argv[1]); /* argument to main  */
        else
            THREADS = 1;                  /* default value     */
        omp_set_num_threads(THREADS);
    #else
        THREADS = 1;                      /* run as sequential */
    #endif


    /* INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION */
    double sigma0 = 1.0;     /* abscissa of convergence */
    unsigned int Nsings = 2; /* number of singularities s_j with Im(s_j) ge 0 */
    /* Singularities */
    double complex SINGS[2] = {-1.0+1.0*I, +1.0+1.0*I};
    /* the following is a gcc extension */
    /* double complex SINGS[2] = {-1.0+1.0i, +1.0+1.0i}; */
    unsigned int MULT[2] = {1,1}; /* polar multiplicities */

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


    /* CALL THE TALBOT SUITE FUNCTION  */
    IFAIL_tot = OMP_Talbot1 (LTfun,sigma0,NTval,t,tol,NUMft,IFAIL,Nsings,SINGS,MULT,Tmin,Tmax);


    /* PRINT RESULTS */
    printf("\nInverting the LT fun at NTval=%d values of t in [%5.2f, %5.2f]\n", NTval,Tmin,Tmax);
    puts("\n\t***   RESULTS OF PARALLEL OMP TALBOT SUITE [function = OMP_Talbot1()]   ***");
    printf("\n\t\tinput tolerance: tol = %5.2e\t\tthreads number: %d\n", tol,THREADS);
    if (NTval < 50)
    {
        if (!IFAIL_tot)
        {
            printf("\n\n       T       F APPROX        IFAIL_tot = %d (no local error)\n", IFAIL_tot);
            for (k=0; k<NTval; k++)
                printf("\n %10.2f   %+e", t[k],NUMft[k]);
        }
        else
        {
            puts("\n\n       T       F APPROX        IFAIL");
            for (k=0; k<NTval; k++)
                printf("\n %10.2f   %+e   %+3d", t[k],NUMft[k],IFAIL[k]);
        }
    }
    puts("\n\n");
    return 0;
}
