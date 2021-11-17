/*
  StencilProbe Heat Equation (Naive version)
  Implements stencil from Chombo's heattut example.
  OpenMP version: Cut over least-stride dimension (nz)
*/
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "stencil.h"


#define NAIVE( LENGTH )                                        \
  double *temp_ptr;                                            \
  int i, j, k, t;                                              \
                                                               \
  for (t = 0; t < timesteps; t++) {                            \
    _Pragma("omp parallel for default(shared) private(i,j,k)") \
    NAIVE_( Anext, A0, LENGTH,                                 \
            LENGTH + PAD(LENGTH), nx - LENGTH - PAD(LENGTH),   \
            LENGTH,               ny - LENGTH,                 \
            LENGTH,               nz - LENGTH )                \
    temp_ptr = A0;                                             \
    A0 = Anext;                                                \
    Anext = temp_ptr;                                          \
  }


/*
  StencilProbe Heat Equation
  Implements 7pt stencil (1)
*/
void StencilProbe_naive_1(double* A0, double* Anext, int nx, int ny, int nz,
                          int tx, int ty, int tz, int timesteps) {
#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  NAIVE( 1 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 13pt stencil (2)
*/
void StencilProbe_naive_2(double* A0, double* Anext, int nx, int ny, int nz,
                          int tx, int ty, int tz, int timesteps) {
#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  NAIVE( 2 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 25pt stencil (4)
*/
void StencilProbe_naive_4(double* A0, double* Anext, int nx, int ny, int nz,
                          int tx, int ty, int tz, int timesteps) {
  NAIVE( 4 );
}


/*
  StencilProbe Heat Equation
  Implements 43pt stencil (7)
*/
void StencilProbe_naive_7(double* A0, double* Anext, int nx, int ny, int nz,
                          int tx, int ty, int tz, int timesteps) {
#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  NAIVE( 7 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 85pt stencil (14)
*/
void StencilProbe_naive_14(double* A0, double* Anext, int nx, int ny, int nz,
                           int tx, int ty, int tz, int timesteps) {
#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  NAIVE( 14 );
#endif
}


#ifdef STENCILTEST
void StencilProbe_naive(double* A0, double* Anext, int nx, int ny, int nz,
                        int tx, int ty, int tz, int timesteps, int length) {
#else
void StencilProbe(double *A0, double *Anext, int nx, int ny, int nz,
                  int tx, int ty, int tz, int timesteps, int length) {
#endif
  switch(length)
  {
    case 1: StencilProbe_naive_1(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 2: StencilProbe_naive_2(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 4: StencilProbe_naive_4(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 7: StencilProbe_naive_7(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 14:StencilProbe_naive_14(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    default: printf("ABORTING: StencilProbe_naive_%d NOT IMPLEMENTED!\n", length); exit(-1);
  }
}

