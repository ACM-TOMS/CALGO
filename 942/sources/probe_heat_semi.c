/*
  StencilProbe Heat Equation (Semi-stencil version)
  Implements stencil from Chombo's heattut example with Semi-stencil.
  OpenMP version: Cut over least-stride dimension (nz), every thread
  computes block sizes of NXxNYxNBZ
*/
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "semi.h"


#ifdef _OPENMP

#define SEMISTENCIL( STENCIL )                                    \
  double *temp_ptr;                                               \
  int i, j, k, t;                                                 \
  int kk, nbz, rnb;                                               \
                                                                  \
  double p2c;                                                     \
  double p2i1, p2i2, p2i3,  p2i4,  p2i5,  p2i6,  p2i7;            \
  double p2i8, p2i9, p2i10, p2i11, p2i12, p2i13, p2i14;           \
  double p2j1, p2j2, p2j3,  p2j4,  p2j5,  p2j6,  p2j7;            \
  double p2j8, p2j9, p2j10, p2j11, p2j12, p2j13, p2j14;           \
  double p2k1, p2k2, p2k3,  p2k4,  p2k5,  p2k6,  p2k7;            \
  double p2k8, p2k9, p2k10, p2k11, p2k12, p2k13, p2k14;           \
                                                                  \
                                                                  \
  _Pragma("omp parallel default(shared) private(kk,nbz,rnb,i,j,k, \
                      t,p2c,p2i1,p2i2,p2i3,p2i4,p2i5,p2i6,p2i7,   \
                      p2i8,p2i9,p2i10,p2i11,p2i12,p2i13,p2i14,    \
                      p2j1,p2j2,p2j3,p2j4,p2j5,p2j6,p2j7,         \
                      p2j8,p2j9,p2j10,p2j11,p2j12,p2j13,p2j14,    \
                      p2k1,p2k2,p2k3,p2k4,p2k5,p2k6,p2k7,         \
                      p2k8,p2k9,p2k10,p2k11,p2k12,p2k13,p2k14)    \
                      firstprivate(A0, Anext, temp_ptr)")         \
  {                                                               \
    nbz = (nz-2*STENCIL) / omp_get_max_threads();                 \
    rnb = (nz-2*STENCIL) % omp_get_max_threads();                 \
    if (rnb > omp_get_thread_num()) {                             \
      nbz++;                                                      \
      kk = STENCIL + nbz*omp_get_thread_num();                    \
    }                                                             \
    else kk = STENCIL + nbz*omp_get_thread_num() + rnb;           \
                                                                  \
    for (t = 0; t < timesteps; t++) {                             \
      SEMISTENCIL_( Anext, A0, STENCIL,                           \
                    STENCIL + PAD(STENCIL), nx - STENCIL - PAD(STENCIL), \
                    STENCIL,                ny - STENCIL,         \
                    kk, MIN(kk+nbz, nz-STENCIL) )                 \
      temp_ptr = A0;                                              \
      A0 = Anext;                                                 \
      Anext = temp_ptr;                                           \
      _Pragma("omp barrier")                                      \
    }                                                             \
  }

#else

#define SEMISTENCIL( STENCIL )                             \
  double *temp_ptr;                                        \
  int i, j, k, t;                                          \
                                                           \
  double p2c;                                              \
  double p2i1, p2i2, p2i3,  p2i4,  p2i5,  p2i6,  p2i7;     \
  double p2i8, p2i9, p2i10, p2i11, p2i12, p2i13, p2i14;    \
  double p2j1, p2j2, p2j3,  p2j4,  p2j5,  p2j6,  p2j7;     \
  double p2j8, p2j9, p2j10, p2j11, p2j12, p2j13, p2j14;    \
  double p2k1, p2k2, p2k3,  p2k4,  p2k5,  p2k6,  p2k7;     \
  double p2k8, p2k9, p2k10, p2k11, p2k12, p2k13, p2k14;    \
                                                           \
                                                           \
  for (t = 0; t < timesteps; t++) {                        \
    SEMISTENCIL_( Anext, A0, STENCIL,                      \
                  STENCIL + PAD(STENCIL), nx - STENCIL - PAD(STENCIL), \
                  STENCIL,                ny - STENCIL,    \
                  STENCIL,                nz - STENCIL )   \
                                                           \
    temp_ptr = A0;                                         \
    A0 = Anext;                                            \
    Anext = temp_ptr;                                      \
  }

#endif


/*
  StencilProbe Heat Equation
  Implements 7pt Semi-stencil (1)
*/
void StencilProbe_semi_1(double* A0, double* Anext, int nx, int ny, int nz,
                         int tx, int ty, int tz, int timesteps) {
#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  SEMISTENCIL( 1 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 13pt Semi-stencil (2)
*/
void StencilProbe_semi_2(double* A0, double* Anext, int nx, int ny, int nz,
                         int tx, int ty, int tz, int timesteps) {
#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  SEMISTENCIL( 2 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 25pt Semi-stencil (4)
*/
void StencilProbe_semi_4(double* A0, double* Anext, int nx, int ny, int nz,
                         int tx, int ty, int tz, int timesteps) {

  SEMISTENCIL( 4 );
}


/*
  StencilProbe Heat Equation
  Implements 43pt Semi-stencil (7)
*/
void StencilProbe_semi_7(double* A0, double* Anext, int nx, int ny, int nz,
                         int tx, int ty, int tz, int timesteps) {
#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  SEMISTENCIL( 7 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 85pt Semi-stencil (14)
*/
void StencilProbe_semi_14(double* A0, double* Anext, int nx, int ny, int nz,
                          int tx, int ty, int tz, int timesteps) {
#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  SEMISTENCIL( 14 );
#endif
}


#ifdef STENCILTEST
void StencilProbe_semi(double* A0, double* Anext, int nx, int ny, int nz,
                       int tx, int ty, int tz, int timesteps, int length) {
#else
void StencilProbe(double *A0, double *Anext, int nx, int ny, int nz,
                  int tx, int ty, int tz, int timesteps, int length) {
#endif
  switch(length)
  {
    case 1: StencilProbe_semi_1(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 2: StencilProbe_semi_2(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 4: StencilProbe_semi_4(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 7: StencilProbe_semi_7(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 14: StencilProbe_semi_14(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    default: printf("ABORTING: StencilProbe_semi_%d NOT IMPLEMENTED!\n", length); exit(-1);
  }
}

