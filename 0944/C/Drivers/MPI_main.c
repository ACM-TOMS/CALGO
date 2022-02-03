/************************        MPI_main.c       *************************
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *             DRIVER PROGRAM FOR THE MPI-BASED IMPLEMENTATION            *
 *                             OF TALBOT SUITE                            *
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

/* Required GNU gcc compiler option:    -std=gnu99 [for clock_gettime()]
   Required         linker   option:    -lrt       [for clock_gettime()]
*/

#include "mpi.h"

#include "../SRC/main_include.h"
#include "../SRC/MPI_Talbot_pack.h" /* for Talbot Suite's MPI implementation */
#include "../SRC/COM_Talbot_pack.h" /* for Talbot Suite's shared functions   */
#include <time.h>                   /* for timing                            */

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


/* PROTOTYPES OF DATA DISTRIBUTION/COLLECTION UTILITY FUNCTIONS */
void data_distribution (unsigned int NTval, double *t, int *NTval_loc, double **t_loc, double **NUMft_loc,
                        int **IFAIL_loc, int **NTvals, int **displs, int Root, MPI_Comm comm);

void data_collection (int IFAIL_tot_loc, int NTval_loc, int *IFAIL_loc, double *NUMft_loc, int *NTvals, int *displs,
                      int *IFAIL_tot, int *IFAIL, double *NUMft, int Root, MPI_Comm comm);

int data_distrib_param (unsigned int NTval, int **NTvals, int **displs, int Root, MPI_Comm comm);

void data_allocate (int NTval_loc, double **t_loc, double **NUMft_loc, int **IFAIL_loc);


int main(int argc, char *argv[])
/* Parameter to the main function:
    1) jFUN   : switch between MPI_Talbot1 (1) and MPI_Talbot2 (2)
*/
{
    int k; /* for-loop index */


    /* SELECTION OF TALBOT SUITE FUNCTION */
    char *Talbot_fun[] = {"MPI_Talbot1",  /* coarse grain parallelism  */
                          "MPI_Talbot2"}; /*  fine  grain parallelism  */
    unsigned int jFUN; /* switch between MPI_Talbot1 (1)
                                     and MPI_Talbot2 (2) */
    if (argc == 2)
        jFUN = (unsigned int)atoi(argv[1]); /* argument to main */
    else
        jFUN = 1;                           /* default value    */


    /* MPI VARIABLES AND INITIALIZATION */
    int myid, nprocs, Root=0;
    MPI_Init (&argc,&argv);
    MPI_Comm_rank (MPI_COMM_WORLD,&myid);
    MPI_Comm_size (MPI_COMM_WORLD,&nprocs);


    /* INPUT DATA ABOUT LAPLACE TRANSFORM FUNCTION */
    char *LTstring  = "F(s) = s/(s^2+9)^2";
    char *ILTstring = "f(t) = t*sin(3*t)/6";
    double sigma0;            /* abscissa of convergence                       */
    unsigned int NsingsTOT;   /* total number of singularities s_j             */
    unsigned int Nsings;      /* number of singularities s_j with Im(s_j) ge 0 */
    double complex SINGS[10]; /* array for singularities                       */
    unsigned int MULT[10];    /* array for polar multiplicities                */
    NsingsTOT = LTsings(&Nsings, SINGS, MULT, &sigma0);

    if (myid == Root)
    {
        printf("\nLaplace Transform function:\t\t%s\tabscissa of convergence: sigma0 = %f\n", LTstring,sigma0);
        printf("\nInv. Laplace Transform function:\t%s\n\n", ILTstring);
        (NsingsTOT == 1  ?  printf("   singularity    and    multiplicity:\n")
                         : printf("   singularities    and    multiplicities:\n"));
        for (k=0; k<NsingsTOT; k++)
            (fabs(cimag(SINGS[k]))>0 ? printf("\ns(%2d) =%+6.2f%+6.2f * I\t\tmult = %d", k+1,creal(SINGS[k]),cimag(SINGS[k]),MULT[k])
                                     : printf("\ns(%2d) =%+6.2f\t\t\tmult = %d", k+1,creal(SINGS[k]),MULT[k]));
        puts("\n");
    }

    /* INPUT DATA ABOUT MULTIPOINT INVERSION PROBLEM              */
    unsigned int NTval=20;     /* number of t values              */
    double Tmin=100, Tmax=500; /* endpoints of the interval for t */
    double tol = 1e-12;        /* tolerance                       */
    double *t;                 /* t array                         */

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


    /* VARIABLES OF THE main FUNCTION */
    struct timespec ts1, ts2; double T0, T1, TOT_TIME, TOT_TIME_loc; /* for timing */
    int NTval_loc, IFAIL_tot_loc, *NTvals, *displs;
    double *t_loc, *NUMft_loc; int *IFAIL_loc;


    /* CALL THE TALBOT SUITE FUNCTION  */
    switch (jFUN)
    {
        case 1 : /* coarse grain parallelism (modified Talbot's method)  */

            /* (1) DISTRIBUTION OF LOCAL DATA (from Root) */
            data_distribution (NTval,t,&NTval_loc,&t_loc,&NUMft_loc,&IFAIL_loc,
                               &NTvals,&displs,Root,MPI_COMM_WORLD);

            clock_gettime(CLOCK_REALTIME, &ts1);

            /* (2) EXECUTE THE LOCAL ALGORITHM     */
            IFAIL_tot_loc = MPI_Talbot1 (LTfun,sigma0,NTval_loc,t_loc,tol,NUMft_loc,IFAIL_loc,
                                         Nsings,SINGS,MULT,Tmin,Tmax);

            /* (3) COMPUTE LOCAL ELAPSED TIMES */
            clock_gettime(CLOCK_REALTIME, &ts2);
            T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec */
            T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
            TOT_TIME_loc = T0 + T1*1.0e-9;

            /* (4) COMBINATION OF LOCAL RESULTS (only in Root) */
            data_collection (IFAIL_tot_loc, NTval_loc, IFAIL_loc, NUMft_loc, NTvals, displs,
                             &IFAIL_tot, IFAIL, NUMft, Root, MPI_COMM_WORLD);
            break;


        case 2 : /*  fine  grain parallelism (classical Talbot's method) */

            /* (1) ALSO NON-Root PROCS ALLOCATE IN/OUT ARRAYS: t, f(t) AND IFAIL */
            if (myid != Root)
                data_allocate (NTval, &t, &NUMft, &IFAIL);

            /* (2) BROADCAST OF t FROM Root */
            MPI_Bcast((void*)t,NTval,MPI_DOUBLE,Root,MPI_COMM_WORLD);

            clock_gettime(CLOCK_REALTIME, &ts1);

            /* (3) EXECUTE THE LOCAL ALGORITHM */
            IFAIL_tot = MPI_Talbot2 (LTfun,sigma0,NTval,t,tol,NUMft,IFAIL,
                                     Nsings,SINGS,MULT,MPI_COMM_WORLD,Root);

            /* (4) COMPUTE LOCAL ELAPSED TIMES */
            clock_gettime(CLOCK_REALTIME, &ts2);
            T0 = (double)(ts2.tv_sec - ts1.tv_sec);   /* sec */
            T1 = (double)(ts2.tv_nsec - ts1.tv_nsec); /* nano sec */
            TOT_TIME_loc = T0 + T1*1.0e-9;
            break;

        default :
            fprintf(stderr, "\n*** No function corresponds to the value %d ***", jFUN);
            exit(1);

    } /* END switch jFUN */


    /* COMPUTE ONLY IN Root THE MAX LOCAL TIME */
    MPI_Reduce(&TOT_TIME_loc,&TOT_TIME,1,MPI_DOUBLE,MPI_MAX,Root,MPI_COMM_WORLD);


    /* PRINT RESULTS ONLY IN Root */
    if (myid == Root)
    {
        double ft, RELERR, ABSERR; char ALFA;
        printf("\nInverting the LT fun at NTval=%d values of t in [%5.2f, %5.2f]\n", NTval,Tmin,Tmax);
        printf("\n\t***   RESULTS OF PARALLEL MPI TALBOT SUITE [function = %s()]   ***\n", Talbot_fun[jFUN-1]);
        printf("\n\t\tinput tolerance: tol = %5.2e\t\tMPI procs number: %d\telapsed time: %e\n", tol,nprocs,TOT_TIME);
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


void data_distribution (unsigned int NTval, double *t, int *NTval_loc, double **t_loc,
                        double **NUMft_loc, int **IFAIL_loc, int **NTvals, int **displs,
                        int Root, MPI_Comm comm)
/*********************************************************************************
    PURPOSE
    =======
    Compute the parameters required for the local data distribution and
    allocate, on all the procs, local arrays for IN/OUT data.

    INPUT PARAMETERS
    ================
    NTval         - Size of the array t.

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
    Collect together local variables in Root.

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
