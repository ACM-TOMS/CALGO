/**     SEQ_Talbot1_SUM_mex.c   -   mex file     **/

#include "mex.h"

#include "../../../TalbotSuiteDE/COM/main_include.h"
#include "../../../TalbotSuiteDE/FUN_DE/SEQ_Talbot_pack_DE.h"



void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/*****      Gateway function      ******
            ================
    IN MATLAB call it as:

        [IFAIL_tot,NUMft,IFAIL] = SEQ_Talbot1_SUM_mex (CONLAM,CONSIG,CONNU,NOPTS,Urows, U, NTval,Tval);
         1         2     3                             1      2      3     4     5      6  7     8

    INPUT PARAMETERS
    ================
    CONLAM,CONSIG,CONNU: geometrical parameters for Talbot's contour.
    NOPTS      : number of columns in U matrix.
    Urows      : number of rows in U matrix.
    U          : row-wise matrix of LT samples of size (Urows, NOPTS).
    NTval      : number of T values.
    Tval       : array of T values.

    OUTPUT PARAMETERS
    =================
    IFAIL_tot  : global error indicator.
    NUMft      : matrix of the numerical Inverse Laplace Transform values. It is
                 vectorized on mesh points: NUMft(h,k) = u(X(h),Y(h),t(k))
    IFAIL      : matrix of the local error indicators. Is is vectorized like NUMft.
**/
{
    unsigned int h, k; /* for loop indices */

    int NARGin = 8, NARGout = 3; /* number of input and output parameters in MATLAB call */

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


    /*  FS is the matrix of LT samples of size (Urows,NOPTS).
        It is a C (row-wise) matrix.
        MATLAB code also returns a row-wise matrix */
    double complex *FS = (double complex *)malloc(Urows*NOPTS*sizeof(double complex));
    if (FS == NULL)
        mexErrMsgTxt("\n***   ERROR in SEQ_Talbot1_SUM_mex: DYNAMIC ALLOCATION OF FS IS FAILED. ***\n");

    for (h=0; h<Urows; h++)     /* row-index */
        for (k=0; k<NOPTS; k++) /* col-index */
            FS[h*NOPTS + k] = *(Ur + h*NOPTS + k) + I* *(Ui + h*NOPTS + k);
            /* both of them are row-wise matrices */


    /* OUTPUT PARAMETERS */
    double *NUMft;         /* numerical f(t) array         */
    int IFAIL_tot, *IFAIL; /* global and local error flags */


    /*  ALLOCATE ALGORITHM'S OUTPUT ARRAYS
        NUMft, IFAIL are C (row-wise) matrices of size (Urows,NTval).
        NUMft(h,k) = u(X(h),Y(h),t(k)),
        where
            (X,Y) contains the vectorized (row-wise) mesh points.
    */
    NUMft = (double *)malloc(Urows*NTval*sizeof(double));
    if (NUMft == NULL)
        mexErrMsgTxt("\n***   ERROR in SEQ_Talbot1_SUM_mex: DYNAMIC ALLOCATION OF NUMft IS FAILED. ***\n");
    IFAIL = (int *)malloc(Urows*NTval*sizeof(int));
    if (IFAIL == NULL)
        mexErrMsgTxt("\n***   ERROR in SEQ_Talbot1_SUM_mex:: DYNAMIC ALLOCATION OF IFAIL IS FAILED. ***\n");


    /* CALL OMP_Talbot1(): coarse grain parallelism (modified Talbot's method) */
    IFAIL_tot = SEQ_TalbotSUM1_DE (CONLAM,CONSIG,CONNU,NOPTS,Urows,FS,NTval,Tval,NUMft,IFAIL);


    /* COPY OUTPUTS TO MATLAB VARIABLES */
    plhs[0] = mxCreateDoubleScalar((double)IFAIL_tot);
    plhs[1] = mxCreateDoubleMatrix(Urows*NTval, 1, mxREAL); /* NUMft */
    plhs[2] = mxCreateDoubleMatrix(Urows*NTval, 1, mxREAL); /* IFAIL */
    double *pNUMft = mxGetPr(plhs[1]),
           *pIFAIL = mxGetPr(plhs[2]);


    /* SINCE MATLAB ALLOCATES A MATRIX BY COLS AND C BY ROWS WE HAVE TO CHANGE
       FROM C TO MATLAB MEMORY ALLOCATION:
                MATLAB        C
                pNUMft(h,k) = NUMft(h,k)
                pIFAIL(h,k) = IFAIL(h,k)
        WHERE
            pNUMft(h,k) = pNUMft[h + k*Urows]   col-wise allocation
            pIFAIL(h,k) = pIFAIL[h + k*Urows]
            NUMft(h,k)  = NUMft[h*NTval + k]    row-wise allocation
            IFAIL(h,k)  = IFAIL[h*NTval + k]
    */
    for (h=0; h<Urows; h++)     /* row-index */
        for (k=0; k<NTval; k++) /* col-index */
        {
            pNUMft[h+k*Urows] = NUMft[h*NTval+k];
            pIFAIL[h+k*Urows] = (double)IFAIL[h*NTval+k];
        }

    free(NUMft); free(IFAIL); free(FS);
    return;
}
