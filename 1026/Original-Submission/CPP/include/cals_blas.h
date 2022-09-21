#ifndef CALS_CALS_BLAS_H
#define CALS_CALS_BLAS_H

#include "omp.h"

#if CALS_MKL

#include "mkl.h"
#define CALS_BACKEND "MKL"

#elif CALS_FLAME

#include "blis.h"
#include "mkl_lapack.h"
#define CALS_BACKEND "BLIS"

#elif CALS_OPENBLAS

#include "cblas.h"

#define CALS_BACKEND "OPENBLAS"

extern "C" {
void dposv_(const char *,
            const long int *,
            const long int *,
            double *,
            const long int *,
            double *,
            const long int *,
            long int *);
void dpotrf_(const char *, const long int *, double *, const long int *, long int *);
};

inline void dposv(const char *uplo,
                  int const *N,
                  int const *Rhs,
                  double *A,
                  int const *Lda,
                  double *B,
                  int const *Ldb,
                  int const *Info) {
  long int const n = *N;
  long int const rhs = *Rhs;
  long int const lda = *Lda;
  long int const ldb = *Ldb;
  long int info = *Info;
  dposv_(uplo, &n, &rhs, A, &lda, B, &ldb, &info);
}

inline void dpotrf(const char *uplo, int const *N, double *A, int const *Lda, int const *Info) {
  long int const n = *N;
  long int const lda = *Lda;
  long int info = *Info;
  dpotrf_(uplo, &n, A, &lda, &info);
}

#elif CALS_MATLAB

#include <cassert>
#include "blas.h"
#include "lapack.h"
#include "mex.h"

#define CALS_BACKEND "MATLAB"

enum CBLAS_ORDER { CblasColMajor = 0, CblasRowMajor = 1 };

enum CBLAS_TRANSPOSE { CblasTrans = 0, CblasNoTrans = 1 };

enum CBLAS_DIAG { CblasNonUnit = 0, CblasUnit = 1 };

enum CBLAS_UPLO { CblasLower = 0, CblasUpper = 1 };

enum CBLAS_SIDE { CblasRight = 0, CblasLeft = 1 };

inline double cblas_dnrm2(long int n_elements, double *data, long int stride) {
  return dnrm2_(&n_elements, data, &stride);
}

inline double cblas_dasum(long int n_elements, double *data, long int stride) {
  return dasum_(&n_elements, data, &stride);
}

inline long int cblas_idamax(long int n_elements, double *data, long int stride) {
  return idamax_(&n_elements, data, &stride) - 1;
}

inline void
cblas_dcopy(long int n_elements, double const *src_data, long int src_stride, double *dst_data, long int dst_stride) {
  dcopy_(&n_elements, src_data, &src_stride, dst_data, &dst_stride);
}

inline void cblas_dscal(long int n_elements, double scalar, double *data, long int stride) {
  dscal_(&n_elements, &scalar, data, &stride);
}

inline void cblas_daxpy(long int n_elements, double scalar, double *x, long int ldx, double *y, long int ldy) {
  daxpy_(&n_elements, &scalar, x, &ldx, y, &ldy);
}

inline void cblas_dgemm(CBLAS_ORDER Order,
                        CBLAS_TRANSPOSE TransA,
                        CBLAS_TRANSPOSE TransB,
                        long int M,
                        long int N,
                        long int K,
                        double alpha,
                        const double *A,
                        long int lda,
                        const double *B,
                        long int ldb,
                        double beta,
                        double *C,
                        long int ldc) {
  assert(Order == CblasColMajor);
  const char *transa = TransA == CblasTrans ? "T" : "N";
  const char *transb = TransB == CblasTrans ? "T" : "N";
  dgemm_(transa, transb, &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

inline void cblas_dgemv(CBLAS_ORDER order,
                        CBLAS_TRANSPOSE TransA,
                        long int M,
                        long int N,
                        double alpha,
                        const double *A,
                        long int lda,
                        const double *X,
                        long int incX,
                        double beta,
                        double *Y,
                        long int incY) {
  assert(order == CblasColMajor);
  const char *transa = TransA == CblasTrans ? "T" : "N";
  dgemv(transa, &M, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
}

inline void cblas_dtrsm(CBLAS_ORDER Order,
                        CBLAS_SIDE Side,
                        CBLAS_UPLO Uplo,
                        CBLAS_TRANSPOSE TransA,
                        CBLAS_DIAG Diag,
                        long int M,
                        long int N,
                        double alpha,
                        const double *A,
                        long int lda,
                        double *B,
                        long int ldb) {
  assert(Order == CblasColMajor);
  const char *transa = TransA == CblasTrans ? "T" : "N";
  const char *side = Side == CblasLeft ? "L" : "R";
  const char *uplo = Uplo == CblasUpper ? "U" : "L";
  const char *diag = Diag == CblasNonUnit ? "N" : "U";
  dtrsm_(side, uplo, transa, diag, &M, &N, &alpha, A, &lda, B, &ldb);
}

inline void dposv(const char *uplo,
                  int const *N,
                  int const *Rhs,
                  double *A,
                  int const *Lda,
                  double *B,
                  int const *Ldb,
                  int const *Info) {
  long int const n = *N;
  long int const rhs = *Rhs;
  long int const lda = *Lda;
  long int const ldb = *Ldb;
  long int info = *Info;
  dposv_(uplo, &n, &rhs, A, &lda, B, &ldb, &info);
}

inline void dpotrf(const char *uplo, int const *N, double *A, int const *Lda, int const *Info) {
  long int const n = *N;
  long int const lda = *Lda;
  long int info = *Info;
  dpotrf_(uplo, &n, A, &lda, &info);
}

#endif

void set_threads(int threads);

int get_threads();

#endif // CALS_CALS_BLAS_H
