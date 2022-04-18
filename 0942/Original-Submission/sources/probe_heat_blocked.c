/*
  StencilProbe Heat Equation (Cache Blocked version, due to Rivera)
  Implements stencil from Chombo's heattut example with cache blocking.
  OpenMP version: each thread computes several TIxTJxNZ blocks
*/
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "stencil.h"


#define TI tx
#ifdef _OPENMP
# define TJ(LENGTH) (((ny - 2*LENGTH)/omp_get_max_threads() < ty) ? (ny - 2*LENGTH)/omp_get_max_threads() : ty)
#else
# define TJ(LENGTH) ty
#endif
#define TK tz


#define RIVERA( LENGTH )                                             \
  /* Fool compiler so it doesn't insert a constant here */           \
  double fac = A0[0];                                                \
  double *temp_ptr;                                                  \
  int i, ii, j, jj, k, t;                                            \
                                                                     \
  ty = TJ(LENGTH);                                                   \
  for (t = 0; t < timesteps; t++) {                                  \
    _Pragma("omp parallel for default(shared) private(jj,ii,i,j,k)") \
    for (jj = LENGTH; jj < ny - LENGTH; jj+=ty) {                    \
      for (ii = LENGTH + PAD(LENGTH); ii < nx - LENGTH - PAD(LENGTH); ii+=TI) { \
        NAIVE_( Anext, A0, LENGTH,                                   \
                ii, MIN(ii+TI, nx-LENGTH-PAD(LENGTH)),               \
                jj, MIN(jj+ty, ny-LENGTH),                           \
                LENGTH, nz - LENGTH )                                \
      }                                                              \
    }                                                                \
    temp_ptr = A0;                                                   \
    A0 = Anext;                                                      \
    Anext = temp_ptr;                                                \
  }


/*
  StencilProbe Heat Equation
  Implements 7pt Blocked stencil (1)
*/
void StencilProbe_rivera_1(double* A0, double* Anext, int nx, int ny, int nz,
                           int tx, int ty, int tz, int timesteps) {
#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  RIVERA( 1 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 13pt Blocked stencil (2)
*/
void StencilProbe_rivera_2(double* A0, double* Anext, int nx, int ny, int nz,
                           int tx, int ty, int tz, int timesteps) {
#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  RIVERA( 2 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 25pt Blocked stencil (4)
*/
void StencilProbe_rivera_4(double* A0, double* Anext, int nx, int ny, int nz,
                           int tx, int ty, int tz, int timesteps) {
  RIVERA( 4 );
}


/*
  StencilProbe Heat Equation
  Implements 43pt Blocked stencil (7)
*/
void StencilProbe_rivera_7(double* A0, double* Anext, int nx, int ny, int nz,
                           int tx, int ty, int tz, int timesteps) {
#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  RIVERA( 7 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 85pt Blocked stencil (14)
*/
void StencilProbe_rivera_14(double* A0, double* Anext, int nx, int ny, int nz,
                            int tx, int ty, int tz, int timesteps) {
#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  RIVERA( 14 );
#endif
}


#ifdef STENCILTEST
void StencilProbe_rivera(double* A0, double* Anext, int nx, int ny, int nz,
                         int tx, int ty, int tz, int timesteps, int length) {
#else
void StencilProbe(double *A0, double *Anext, int nx, int ny, int nz,
                  int tx, int ty, int tz, int timesteps, int length) {
#endif
  switch(length)
  {
    case 1: StencilProbe_rivera_1(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 2: StencilProbe_rivera_2(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 4: StencilProbe_rivera_4(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 7: StencilProbe_rivera_7(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 14:StencilProbe_rivera_14(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    default: printf("ABORTING: StencilProbe_rivera_%d NOT IMPLEMENTED!\n", length); exit(-1);
  }
}

