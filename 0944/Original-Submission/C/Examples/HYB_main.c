/************************        HYB_main.c       *************************
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *                 SAMPLE MAIN for HYB_Talbot3() function                 *
 *        IN THE MPI/OMP-BASED HYB IMPLEMENTATION OF TALBOT SUITE         *
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *  >>>>>>>>>>>>>      VERSION 1.0    Sept 13th, 2012         <<<<<<<<<<  *
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

#include "mpi.h"

#include "../SRC/main_include.h"
#include "../SRC/HYB_Talbot_pack.h" /* for Talbot Suite's HYB implementation */
#include "../SRC/COM_Talbot_pack.h" /* for Talbot Suite's shared functions   */

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


/* PROTOTYPES OF DATA DISTRIBUTION UTILITY FUNCTIONS */
void data_distribution (unsigned int NTval, double *t, int *NTval_loc, double **t_loc, double **NUMft_loc,
                        int **IFAIL_loc, int **NTvals, int **displs, int Root, MPI_Comm comm);

void data_collection (int IFAIL_tot_loc, int NTval_loc, int *IFAIL_loc, double *NUMft_loc, int *NTvals, int *displs,
                      int *IFAIL_tot, int *IFAIL, double *NUMft, int Root, MPI_Comm comm);

void data_allocate (int NTval_loc, double **t_loc, double **NUMft_loc, int **IFAIL_loc);

int data_distrib_param (unsigned int NTval, int **NTvals, int **displs, int Root, MPI_Comm comm);


int main(int argc, char *argv[])
{
    int k; /* for-loop index */


    /* MPI VARIABLES AND INITIALIZATION */
    int nprocs, myid, Root=0;
    MPI_Init (&argc,&argv);
    MPI_Comm_rank (MPI_COMM_WORLD,&myid);
    MPI_Comm_size (MPI_COMM_WORLD,&nprocs);


    /* OpenMP ENVIRONMENT SETTINGS */
    int THREADS; /* number of threads in parallel sections */
    #ifdef _OPENMP
        omp_set_dynamic(0); /* static adjustment of thread number */
        /* SET THE NUMBER OF THREADS */
        if (argc == 2)
            THREADS = (int)atoi(argv[1]); /* argument to main */
        else
            THREADS = 1;                  /* default value    */
        omp_set_num_threads(THREADS);
    #else
        THREADS = 1;                      /* run as sequential */
    #endif


    /* INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION (on all the procs) */
    double sigma0 = 1.0;     /* abscissa of convergence */
    unsigned int Nsings = 2; /* number of singularities s_j with Im(s_j) ge 0 */
    /* Singularities */
    double complex SINGS[2] = {-1.0+1.0*I, +1.0+1.0*I};
    /* the following is a gcc extension */
    /* double complex SINGS[2] = {-1.0+1.0i, +1.0+1.0i}; */
    unsigned int MULT[2] = {1,1}; /* polar multiplicities */


    /* INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM (on all the procs) */
    unsigned int NTval=20;     /* number of t values                    */
    double Tmin=100, Tmax=500; /* endpoints of the interval for t       */
    double tol = 1e-12;        /* tolerance                             */
    double *t;                 /* t array                               */


    /* OUTPUT DATA */
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot, *IFAIL; /* global and local error flags */

    if (myid == Root)
    {
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
            *(t+k)  = Tmin + k*t_step;

        /* ALLOCATE ALGORITHM'S OUTPUT ARRAYS */
        NUMft = (double *)malloc(NTval*sizeof(double));
        if (NUMft == NULL)
            {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n"); exit(1);}

        IFAIL = (int *)malloc(NTval*sizeof(int));
        if (IFAIL == NULL)
            {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n"); exit(1);}
    }


    /* (1) DISTRIBUTION OF LOCAL DATA (from Root) */
    int NTval_loc, IFAIL_tot_loc, *NTvals, *displs;
    double *t_loc, *NUMft_loc; int *IFAIL_loc;

    data_distribution (NTval,t,&NTval_loc,&t_loc,&NUMft_loc,&IFAIL_loc,&NTvals,&displs,Root,MPI_COMM_WORLD);

    /* (2) EXECUTE THE LOCAL ALGORITHM */
    /*     CALL THE PARALLEL HYB FUNCTION - modified Talbot's method
           (MPI coarse grain parallelism / OMP fine grain parallelism) */
    IFAIL_tot_loc = HYB_Talbot3 (LTfun,sigma0,NTval_loc,t_loc,tol,NUMft_loc,IFAIL_loc,Nsings,SINGS,MULT,Tmin,Tmax);

    /* (3) COMBINATION OF LOCAL RESULTS (only in Root) */
    data_collection (IFAIL_tot_loc, NTval_loc, IFAIL_loc, NUMft_loc, NTvals, displs,
                             &IFAIL_tot, IFAIL, NUMft, Root, MPI_COMM_WORLD);


    /* PRINT RESULTS ONLY IN Root */
    if (myid == Root)
    {
        printf("\nInverting the LT fun at NTval=%d values of t in [%5.2f, %5.2f]\n", NTval,Tmin,Tmax);
        puts("\n\t***   RESULTS OF PARALLEL HYB TALBOT SUITE [function = HYB_Talbot3()]   ***");
        printf("\n\tinput tolerance: tol = %5.2e\t\tMPI procs: %d\t\tOMP threads: %d\n", tol,nprocs,THREADS);
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
    }
    MPI_Finalize();
    return 0;
}


void data_allocate (int NTval_loc, double **t_loc, double **NUMft_loc, int **IFAIL_loc)
/*********************************************************************************
    PURPOSE
    =======
    Allocate, on all the procs, local arrays for IN/OUT data.

    INPUT PARAMETERS
    ================
    NTval_loc     - Size of local arrays to be allocated.

    OUTPUT PARAMETERS
    =================
    t_loc, NUMft_loc, IFAIL_loc - Local arrays allocated by the function.
 *********************************************************************************/
{
    *t_loc = (double *)malloc(NTval_loc*sizeof(double));
    if (*t_loc == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF t_loc IS FAILED. ***\n"); exit(1);}

    *NUMft_loc = (double *)malloc(NTval_loc*sizeof(double));
    if (*NUMft_loc == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF NUMft_loc IS FAILED. ***\n"); exit(1);}

    *IFAIL_loc = (int *)malloc(NTval_loc*sizeof(int));
    if (*IFAIL_loc == NULL)
        {fprintf(stderr, "\n***   ERROR: DYNAMIC ALLOCATION OF IFAIL_loc IS FAILED. ***\n"); exit(1);}
}


int data_distrib_param (unsigned int NTval, int **NTvals, int **displs, int Root, MPI_Comm comm)
/*********************************************************************************
    PURPOSE
    =======
    Compute the parameters required for the local data distribution.

    INPUT PARAMETERS
    ================
    NTval         - Size of data arrays depending on the t values.

    OUTPUT PARAMETERS
    =================
    NTval_loc      - (int) local size of data arrays.
                    This value is returned by the function itself in all the procs.

    NTvals, displs - (int arrays) NTvals is the array containing the local sizes
                    NTval_loc for further data distribution; displs is the array
                    containing the displacements required by MPI_Scatterv to distribute
                    t values locally. They are returned only in Root.
 *********************************************************************************/
{
    int nprocs, myid, rem, k;
    int NTval_loc, *p_NTvals=NULL;

    MPI_Comm_size (comm,&nprocs);
    MPI_Comm_rank (comm,&myid);

    /* All procs compute the sizes of local t values */
    NTval_loc = NTval/nprocs;   rem = NTval%nprocs;
    if (myid < rem)
        (NTval_loc)++;

    /* Root receives in NTvals the sizes (NTval_loc) of local arrays */
    if (myid == Root)
    {
        *NTvals = (int*)calloc(nprocs, sizeof(int));
        p_NTvals = *NTvals;
    }
    MPI_Gather((void*)&NTval_loc,1,MPI_UNSIGNED,(void*)p_NTvals,1,MPI_UNSIGNED,Root,comm);

    /* Root computes the start address (displacement) of local array in the global array of t values */
    if(myid == Root)
    {
        *displs = (int*)malloc(nprocs*sizeof(int));
        *(*displs) = 0;
        for (k=1; k<nprocs; k++)
            *(*displs+k) = *(*displs+k-1) + *(*NTvals+k-1);
    }

    return NTval_loc;
}




void data_distribution (unsigned int NTval, double *t, int *NTval_loc, double **t_loc, double **NUMft_loc, int **IFAIL_loc, int **NTvals, int **displs, int Root, MPI_Comm comm)
/*********************************************************************************
    PURPOSE
    =======
    Compute the parameters required for the local data distribution and
    allocate, on all the procs, local arrays for IN/OUT data.

    INPUT PARAMETERS
    ================
    NTval         - Size of data arrays depending on the t values.

    t             - Array of t values.

    comm          - (MPI communicator) identifies the group of processes
                     involved in the process.

    Root          - (integer) rank of Root process.


    OUTPUT PARAMETERS
    =================
    NTval_loc      - (int) local size of data arrays.

    NTvals, displs - (int arrays) NTvals is the array containing the local sizes
                    NTval_loc for data distribution; displs is the array containing
                    the displacements required by MPI_Scatterv to distribute
                    t values locally.

    t_loc, NUMft_loc, IFAIL_loc - Local arrays allocated by the function.

 *********************************************************************************/
{
    /* DISTRIBUTION OF LOCAL DATA (from Root) */
    *NTval_loc = data_distrib_param (NTval, NTvals, displs, Root, comm);
    data_allocate (*NTval_loc, t_loc, NUMft_loc, IFAIL_loc);
    MPI_Scatterv((void*)t, *NTvals, *displs, MPI_DOUBLE, (void*)(*t_loc), *NTval_loc, MPI_DOUBLE, Root, comm);
}


void data_collection (int IFAIL_tot_loc, int NTval_loc, int *IFAIL_loc, double *NUMft_loc, int *NTvals, int *displs,
                      int *IFAIL_tot, int *IFAIL, double *NUMft, int Root, MPI_Comm comm)
/*********************************************************************************
    PURPOSE
    =======
    Collect together local variables.

    INPUT PARAMETERS
    ================
    IFAIL_tot_loc  - (integer) local global error indicator.

    NTval_loc      - (integer) local size of data arrays (in all the procs).

    IFAIL_loc      - (integer array of size NTval_loc) local error indicators.

    NUMft_loc      - (double array of size NTval_loc) local f(t) values.

    NTvals, displs - (int arrays) NTvals is the array containing the local sizes
                     NTval_loc for data collection; displs is the array containing
                     the displacements required by MPI_Gatherv to collect local
                     variables in Root.

    comm           - (MPI communicator) identifies the group of processes
                     involved in the process.

    Root           - (integer) rank of Root process.


    OUTPUT PARAMETERS
    =================
    IFAIL_tot      - (integer) global error indicator computed as a logical or
                     among the values of IFAIL_tot_loc.

    IFAIL          - (integer array) array gathered to Root process from local IFAIL_loc.

    NUMft          - (double array) array gathered to Root process from local NUMft_loc.

 *********************************************************************************/
{
    MPI_Reduce(&IFAIL_tot_loc, IFAIL_tot, 1, MPI_INT, MPI_LOR, Root, comm);
    MPI_Gatherv((void*)IFAIL_loc, NTval_loc, MPI_INT,(void*)IFAIL, NTvals, displs, MPI_INT, Root, comm);
    MPI_Gatherv((void*)NUMft_loc, NTval_loc, MPI_DOUBLE,(void*)NUMft, NTvals, displs, MPI_DOUBLE, Root, comm);
}
