/*  sdpa_algebra.h

LAPACK+BLAS definitions wrapper

Define macros to mangle the given C identifier (in lower and upper
case), which must not contain underscores, for linking with Fortran.

*/

#ifndef __sdpa_algebra_h__
#define __sdpa_algebra_h__

#define F77_RET_I int
#define F77_RET_D double

#if defined(__APPLE__) // Dirty...
#define F77_FUNC(name,NAME) name ## _
#endif

#define dtrsm_f77  F77_FUNC (dtrsm,  DTRSM)
#define dsyrk_f77  F77_FUNC (dsyrk,  DSYRK)
#define dcopy_f77  F77_FUNC (dcopy,  DCOPY)
#define daxpy_f77  F77_FUNC (daxpy,  DAXPY)
#define dgemm_f77  F77_FUNC (dgemm,  DGEMM)
#define dgemv_f77  F77_FUNC (dgemv,  DGEMV)
#define dscal_f77  F77_FUNC (dscal,  DSCAL)
#define dtrsv_f77  F77_FUNC (dtrsv,  DTRSV)
#define dtrmv_f77  F77_FUNC (dtrmv,  DTRMV)
#define ddot_f77   F77_FUNC (ddot,   DDOT)
#define dtrmm_f77  F77_FUNC (dtrmm,  DTRMM)
#define ilaenv_f77 F77_FUNC (ilaenv, ILAENV)
#define dsteqr_f77 F77_FUNC (dsteqr, DSTEQR)
#define dsyev_f77  F77_FUNC (dsyev,  DSYEV)
#define dpotrf_f77 F77_FUNC (dpotrf, DPORTRF)

extern "C"
{
// BLAS
  F77_RET_I  dtrsm_f77
      (char* side, char* uplo, char* trans, char* diag,
       int* M, int* N,
       double* alpha,
       double* A, int* lda,
       double* B, int* ldb, int side_len,
       int uplo_len, int trans_len, int diag_len);

  F77_RET_I  dsyrk_f77
      (char* uplo, char* trans, int* N, int* K,
       double* alpha,
       double* A, int* lda,
       double* beta,
       double* C, int* ldc, int uplo_len, int trans_len);

  F77_RET_I  dcopy_f77
      (int* N,
       double* X, int* incX,
       double* Y, int* incY);

  F77_RET_I  daxpy_f77
      (int* N,
       double* alpha,
       double* X, int* incX,
       double* Y, int* incY);

  F77_RET_I  dgemm_f77
      (char* transA, char* transB, int* M, int* N, int* K,
       double* alpha,
       double* A, int* lda,
       double* B, int* ldb,
       double* beta,
       double* C, int* ldc, int transA_len, int transB_len);

  F77_RET_I  dgemv_f77
      (char* trans, int* M, int* N,
       double* alpha,
       double* A, int* lda,
       double* X, int* incX,
       double* beta,
       double* Y, int* incY, int trans_len);

  F77_RET_I  dscal_f77
      (int* N,
       double* alpha,
       double* X, int* incX);

  F77_RET_I  dtrsv_f77
      (char* uplo, char* trans, char* diag, int* N,
       double* A, int* lda,
       double* X, int* incX, int uplo_len,
       int trans_len, int diag_len);

  F77_RET_I  dtrmv_f77
      (char* uplo, char *trans, char* diag, int *N,  
       double *A, int *lda, 
       double *X, int *incX, int uplo_len, int trans_len, int diag_len);

  F77_RET_D  ddot_f77
      (int* N, double* X, int* incX, double* Y, int* incY);

  F77_RET_I  dtrmm_f77
      (char* side, char* uplo, char* trans, char* diag, 
       int* M, int* N,
       double* alpha,
       double* A, int* lda,
       double* B, int* ldb, int side_len, int uplo_len,
       int trans_len, int diag_len);

// LAPACK

  F77_RET_I  ilaenv_f77
      (int *ispec, char *name, char *opts, int *n1, 
	int *n2, int *n3, int *n4, int name_len, int opts_len);

  F77_RET_I  dsteqr_f77
      (char *compz, int *n, double *d, 
	double *e, double *z, int *ldz, double *work, 
	int *info, int compz_len);

  F77_RET_I  dsyev_f77
      (char *jobz, char *uplo, int *n, double *a,
        int *lda, double *w, double *work, int *lwork, 
	int *info, int jobz_len, int uplo_len);

  F77_RET_I  dpotrf_f77
     (char *uplo, int *n, double *a, int *lda,
      int *info, int uplo_len);

}

#endif // __sdpa_algebra_h__
