/**     OMP_Talbot11_SUM_mex.c   -   mex file     **/

#include "mex.h"

#include "../../../TalbotSuiteDE/COM/main_include.h"
#include "../../../TalbotSuiteDE/FUN_DE/OMP_Talbot_pack_DE.h"

#ifdef _OPENMP
    #include <omp.h>
#endif



void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/*****      Gateway function      ******
            ================
    IN MATLAB call it as:

        [IFAIL_tot,NUMft,IFAIL] = OMP_Talbot1_SUM_mex (CONLAM,CONSIG,CONNU,NOPTS,Urows, U, NTval,Tval,THREADS);
         1         2     3                             1      2      3     4     5      6  7     8    9

    INPUT PARAMETERS
    ================
    CONLAM,CONSIG,CONNU: geometrical parameters for Talbot's contour.
    NOPTS      : number of columns in U matrix.
    Urows      : number of rows in U matrix.
    U          : row-wise matrix of LT samples of size (Urows, NOPTS).
    NTval      : number of T values.
    Tval       : array of T values.
    THREADS    : number of parallel threads

    OUTPUT PARAMETERS
    =================
    IFAIL_tot  : global error indicator.
    NUMft      : matrix of the numerical Inverse Laplace Transform values. It is
                 vectorized on mesh points: NUMft(h,k) = u(X(h),Y(h),t(k))
    IFAIL      : matrix of the local error indicators. Is is vectorized like NUMft.
**/
{
    unsigned int h, j, k; /* for loop indices */

    int NARGin = 9, NARGout = 3; /* number of input and output parameters in MATLAB call */

    /* Check for the proper number of input/output arguments. */
    if (nrhs != NARGin)
        mexErrMsgTxt("Wrong number of input parameters.");

    if (nlhs != NARGout)
        mexErrMsgTxt("Wrong number of output parameters.");


    /* EXTRACT INPUT PARAMETERS  */
    double CONLAM = mxGetScalar(prhs[0]),
           CONSIG = mxGetScalar(prhs[1]),
           CONNU  = mxGetScalar(prhs[2]);
    unsigned int NOPTS = (unsigned int)mxGetScalar(prhs[3]),
                 Urows = (unsigned int)mxGetScalar(prhs[4]);
    unsigned int NTval = (unsigned int)mxGetScalar(prhs[6]); /* number of t values */
    double *Tval = mxGetPr(prhs[7]);                         /* t array */
    int        THREADS = (int)mxGetScalar(prhs[8]);          /* OpenMP parallel threads */

    /*  EXTRACT THE 6th PARAMETER: U = Ur + i*Ui
        U is the row-wise matrix of LT samples of size (Urows,NOPTS).
                U(h,k) = U(X(h),Y(h),S(k))
        where
            (X,Y) contains the vectorized (row-wise) mesh points.
            S contains the Talbot contour points.
     */
    double  *Ur, *Ui;
    Ur = mxGetPr(prhs[5]);
    Ui = mxGetPi(prhs[5]);


    unsigned int NTOT = Urows*NOPTS;

    /*  FS is the matrix of LT samples of size (Urows,NOPTS): FS(h,k) = U(X(h),Y(h),S(k)).
        It is a C (row-wise) matrix.  MATLAB code also provides a row-wise matrix U */
    double complex *FS = (double complex *)malloc(NTOT*sizeof(double complex));
    if (FS == NULL)
        mexErrMsgTxt("\n***   ERROR in SEQ_Talbot1_SUM_mex: DYNAMIC ALLOCATION OF FS IS FAILED. ***\n");
    #pragma omp parallel for    default   (shared)   \
                                private   (h)        \
                                num_threads (THREADS)
    for (h=0; h<NTOT; h++)
        FS[h] = Ur[h] + I*Ui[h];
        /* both of them are row-wise matrices */


    /* OUTPUT PARAMETERS */
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot, *IFAIL; /* global and local error flags */


    /*  ALLOCATE ALGORITHM'S OUTPUT ARRAYS
        NUMft, IFAIL are C (row-wise) matrices of size (Urows,NTval).
        NUMft(h,k) = u(X(h),Y(h),t(k)),
        where (X,Y) contain the vectorized (row-wise) mesh points. */
    NUMft = (double *)malloc(Urows*NTval*sizeof(double));
    if (NUMft == NULL)
        mexErrMsgTxt("\n***   ERROR in SEQ_Talbot1_SUM_mex: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n");
    IFAIL = (int *)malloc(Urows*NTval*sizeof(int));
    if (IFAIL == NULL)
        mexErrMsgTxt("\n***   ERROR in SEQ_Talbot1_SUM_mex:: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n");


    /* CALL OMP_Talbot1(): coarse grain parallelism (modified Talbot's method) */
    IFAIL_tot = OMP_TalbotSUM11_DE (CONLAM,CONSIG,CONNU,NOPTS,Urows,FS,NTval,Tval,NUMft,IFAIL,THREADS);


    /* COPY OUTPUTS TO MATLAB VARIABLES */
    plhs[0] = mxCreateDoubleScalar((double)IFAIL_tot);
    plhs[1] = mxCreateDoubleMatrix(Urows*NTval, 1, mxREAL); /* NUMft */
    plhs[2] = mxCreateDoubleMatrix(Urows*NTval, 1, mxREAL); /* IFAIL */
    double *pNUMft = mxGetPr(plhs[1]),
           *pIFAIL = mxGetPr(plhs[2]);


    /* SINCE MATLAB ALLOCATES A MATRIX BY COLS AND C BY ROWS WE HAVE TO CHANGE
       FROM C TO MATLAB MEMORY ALLOCATION:
                MATLAB        C
                pNUMft(h,k) = NUMft(h,k)    h=0, ..., Nrows-1
                pIFAIL(h,k) = IFAIL(h,k)    k=0, ..., NTval-1
        WHERE
                pNUMft(h,k) = pNUMft[h + k*Urows]   col-wise allocation
                pIFAIL(h,k) = pIFAIL[h + k*Urows]
                NUMft(h,k)  = NUMft[h*NTval + k]    row-wise allocation
                IFAIL(h,k)  = IFAIL[h*NTval + k]

        h and k MAY BE COMPUTED AS:
                h = j%Nrows   ==>   k = j/Nrows
        OR
                k = j%NTval   ==>   h = j/NTval
        FOR j=0, ..., NTOT-1
    */
    NTOT = Urows*NTval;
    unsigned int NTOTloc, STARTloc, ENDloc;
    int          mod, myid;

    #pragma omp parallel for    default   (shared)   \
                                private   (h,j,k)    \
                                num_threads (THREADS)
    for (j=0; j<NTOT; j++)
    {
        h = j%Urows;
        k = j/Urows;
        pNUMft[h+k*Urows] = NUMft[h*NTval+k];
        pIFAIL[h+k*Urows] = (double)IFAIL[h*NTval+k];
    }

    free(NUMft); free(IFAIL); free(FS);
    return;
}
