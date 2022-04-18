/*  Time skewing stencil code
 *  Kaushik Datta (kdatta@cs.berkeley.edu)
 *  University of California Berkeley
 *
 *  This code implements the time skewing method.  The cache blocks need to be
 *  traversed in a specific order for the algorithm to work properly.
 *
 *  NOTE: The number of iterations can only be up to one greater than the
 *  smallest cache block dimension.  If you wish to do more iterations, there
 *  are two options:
 *    1.  Make the smallest cache block dimension larger.
 *    2.  Split the number of iterations into smaller runs where each run
 *        conforms to the above rule.
 */
/*
  StencilProbe Heat Equation (Timeskew Semi-stencil version)
  Implements stencil from Chombo's heattut example with timeskew & Semi-stencil.
  Modified by Raúl de la Cruz (delacruz@bsc.es)
  Barcelona Supercomputing Center

  Based on Timeskewing algorithm from "Time Skewing: A Value-Based Approach to
  Optimizing for Memory Locality" of John McCalpin and David Wonnacott.
  OpenMP version: based on "Time Skewing for Parallel Computers" paper of David
  Wonnacott. Parallelism is created in space cuts of Z dimension (every thread
  computes blocks of TXxTYxTZ size), and three types of parallelepipeds (light,
  white and black) are computed by each thread to preserve dependencies.
*/
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "semi.h"


#ifdef _OPENMP

#define TIMESKEW_SEMI( LENGTH )                                      \
  double p2c;                                                        \
  double p2i1, p2i2, p2i3,  p2i4,  p2i5,  p2i6,  p2i7;               \
  double p2i8, p2i9, p2i10, p2i11, p2i12, p2i13, p2i14;              \
  double p2j1, p2j2, p2j3,  p2j4,  p2j5,  p2j6,  p2j7;               \
  double p2j8, p2j9, p2j10, p2j11, p2j12, p2j13, p2j14;              \
  double p2k1, p2k2, p2k3,  p2k4,  p2k5,  p2k6,  p2k7;               \
  double p2k8, p2k9, p2k10, p2k11, p2k12, p2k13, p2k14;              \
                                                                     \
  double fac = A0[0];                                                \
  double *temp_ptr;                                                  \
  double *myA0, *myAnext;                                            \
                                                                     \
  int neg_x_slope, neg_y_slope, neg_z_slope;                         \
  int pos_x_slope, pos_y_slope, pos_z_slope;                         \
  int blockMin_x, blockMin_y, blockMin_z;                            \
  int blockMax_x, blockMax_y, blockMax_z;                            \
  int ii, jj, kk, i, j, k, t;                                        \
  int txx, tyy, tzz;                                                 \
  int ttz, tk, rk, ik;                                               \
  int zz, base;                                                      \
                                                                     \
                                                                     \
  /* Base size for parallelepipeds */                                \
  base = LENGTH*(timesteps-1)*2+1;                                   \
                                                                     \
  /* Recompute tz block size to be tz*#threads if not enough    */   \
  /* work for each thread. Also tz >= than base (critical size) */   \
  tk = (nz - 2*LENGTH)/omp_get_max_threads();                        \
  rk = (nz - 2*LENGTH)%omp_get_max_threads();                        \
  if (tk < tz) ttz = nz - 2*LENGTH;                                  \
  else {                                                             \
    tk = tz; rk = 0;                                                 \
    ttz = tz * omp_get_max_threads();                                \
  }                                                                  \
  if (base > tk) {                                                   \
    printf("Error: base (%d) must be <= than tz (%d)\n", base, tk ); \
    return;                                                          \
  }                                                                  \
                                                                     \
  for (kk= LENGTH; kk < nz-LENGTH; kk+= ttz) {                       \
    neg_z_slope = LENGTH;                                            \
    pos_z_slope = -LENGTH;                                           \
                                                                     \
    for (jj= LENGTH; jj < ny-LENGTH; jj+= ty) {                      \
      neg_y_slope = LENGTH;                                          \
      pos_y_slope = -LENGTH;                                         \
      tyy = MIN(ty, ny-jj-LENGTH);                                   \
                                                                     \
      if (jj == LENGTH) {                                            \
        neg_y_slope = 0;                                             \
      }                                                              \
      if (jj == ny-tyy-LENGTH) {                                     \
        pos_y_slope = 0;                                             \
      }                                                              \
      for (ii= LENGTH + PAD(LENGTH); ii < nx-LENGTH - PAD(LENGTH); ii+= tx) { \
        neg_x_slope = LENGTH;                                        \
        pos_x_slope = -LENGTH;                                       \
        txx = MIN(tx, nx-ii-LENGTH - PAD(LENGTH));                   \
                                                                     \
        if (ii == LENGTH + PAD(LENGTH)) {                            \
          neg_x_slope = 0;                                           \
        }                                                            \
        if (ii == nx-txx-LENGTH - PAD(LENGTH)) {                     \
          pos_x_slope = 0;                                           \
        }                                                            \
                                                                     \
        _Pragma("omp parallel default(shared) private(i,j,k,t,zz,tzz,\
                            blockMin_x,blockMin_y,blockMin_z,tz,     \
                            blockMax_x,blockMax_y,blockMax_z,ik,     \
                            p2c,p2i1,p2i2,p2i3,p2i4,p2i5,p2i6,p2i7,  \
                            p2i8,p2i9,p2i10,p2i11,p2i12,p2i13,p2i14, \
                            p2j1,p2j2,p2j3,p2j4,p2j5,p2j6,p2j7,      \
                            p2j8,p2j9,p2j10,p2j11,p2j12,p2j13,p2j14, \
                            p2k1,p2k2,p2k3,p2k4,p2k5,p2k6,p2k7,      \
                            p2k8,p2k9,p2k10,p2k11,p2k12,p2k13,p2k14) \
                          firstprivate(myA0,myAnext,temp_ptr,        \
                                       neg_z_slope,pos_z_slope)")    \
        {                                                            \
          /* Ensure good load-balance among threads */               \
          tz = tk + ((rk > omp_get_thread_num()) ? 1 : 0);           \
          ik = tk * omp_get_thread_num() +                           \
                ((rk > omp_get_thread_num()) ? omp_get_thread_num() : rk); \
                                                                     \
          zz  = kk + ik;                                             \
          tzz = MIN(tz, nz-zz-LENGTH);                               \
                                                                     \
          /* Don't compute if tz block is smaller than base */       \
          /* Previous domain will compute its domain part.  */       \
          if (tzz < base) tzz = 0;                                   \
          else if (nz-LENGTH-zz-tzz < base) tzz += nz-LENGTH-zz-tzz; \
                                                                     \
          if (zz == LENGTH) {                                        \
            neg_z_slope = 0;                                         \
          }                                                          \
          if (zz == nz-tzz-LENGTH) {                                 \
            pos_z_slope = 0;                                         \
          }                                                          \
                                                                     \
          myA0 = A0;                                                 \
          myAnext = Anext;                                           \
                                                                     \
          /* Compute light parallelepiped */                         \
          for (t= 0; t < timesteps && tzz > 0; t++) {                \
            blockMin_x = MAX(LENGTH + PAD(LENGTH), ii - t * neg_x_slope); \
            blockMin_y = MAX(LENGTH, jj - t * neg_y_slope);          \
            blockMin_z = MAX(LENGTH, zz + tzz - base + t * LENGTH);  \
                                                                     \
            blockMax_x = MAX(LENGTH + PAD(LENGTH), ii + txx + t * pos_x_slope); \
            blockMax_y = MAX(LENGTH, jj + tyy + t * pos_y_slope);    \
            blockMax_z = MIN(nz-LENGTH, zz + tzz + t * pos_z_slope); \
                                                                     \
            SEMISTENCIL_( myAnext, myA0, LENGTH,                     \
                          blockMin_x, blockMax_x,                    \
                          blockMin_y, blockMax_y,                    \
                          blockMin_z, blockMax_z )                   \
                                                                     \
            temp_ptr = myA0;                                         \
            myA0 = myAnext;                                          \
            myAnext = temp_ptr;                                      \
          }                                                          \
                                                                     \
          myA0 = A0;                                                 \
          myAnext = Anext;                                           \
                                                                     \
          /* Compute white parallelepiped */                         \
          for (t= 0; t < timesteps && tzz > 0; t++) {                \
            blockMin_x = MAX(LENGTH + PAD(LENGTH), ii - t * neg_x_slope); \
            blockMin_y = MAX(LENGTH, jj - t * neg_y_slope);          \
            blockMin_z = MAX(LENGTH, zz + t * LENGTH);               \
                                                                     \
            blockMax_x = MAX(LENGTH + PAD(LENGTH), ii + txx + t * pos_x_slope); \
            blockMax_y = MAX(LENGTH, jj + tyy + t * pos_y_slope);    \
            blockMax_z = MAX(LENGTH, zz + tzz - base + t * LENGTH);  \
                                                                     \
            SEMISTENCIL_( myAnext, myA0, LENGTH,                     \
                          blockMin_x, blockMax_x,                    \
                          blockMin_y, blockMax_y,                    \
                          blockMin_z, blockMax_z )                   \
                                                                     \
            temp_ptr = myA0;                                         \
            myA0 = myAnext;                                          \
            myAnext = temp_ptr;                                      \
          }                                                          \
                                                                     \
          myA0 = A0;                                                 \
          myAnext = Anext;                                           \
                                                                     \
          /* Syncrhonize threads before dark region computation */   \
          _Pragma("omp barrier")                                     \
                                                                     \
          /* Compute dark parallelepiped */                          \
          for (t= 0; t < timesteps && tzz > 0; t++) {                \
            blockMin_x = MAX(LENGTH + PAD(LENGTH), ii - t * neg_x_slope); \
            blockMin_y = MAX(LENGTH, jj - t * neg_y_slope);          \
            blockMin_z = MAX(LENGTH, zz - t * neg_z_slope);          \
                                                                     \
            blockMax_x = MAX(LENGTH + PAD(LENGTH), ii + txx + t * pos_x_slope); \
            blockMax_y = MAX(LENGTH, jj + tyy + t * pos_y_slope);    \
            blockMax_z = MAX(LENGTH, zz + t * LENGTH);               \
                                                                     \
            SEMISTENCIL_( myAnext, myA0, LENGTH,                     \
                          blockMin_x, blockMax_x,                    \
                          blockMin_y, blockMax_y,                    \
                          blockMin_z, blockMax_z )                   \
                                                                     \
            temp_ptr = myA0;                                         \
            myA0 = myAnext;                                          \
            myAnext = temp_ptr;                                      \
          }                                                          \
        }                                                            \
      }                                                              \
    }                                                                \
  }

#else

#define TIMESKEW_SEMI( LENGTH )                                      \
  double p2c;                                                        \
  double p2i1, p2i2, p2i3,  p2i4,  p2i5,  p2i6,  p2i7;               \
  double p2i8, p2i9, p2i10, p2i11, p2i12, p2i13, p2i14;              \
  double p2j1, p2j2, p2j3,  p2j4,  p2j5,  p2j6,  p2j7;               \
  double p2j8, p2j9, p2j10, p2j11, p2j12, p2j13, p2j14;              \
  double p2k1, p2k2, p2k3,  p2k4,  p2k5,  p2k6,  p2k7;               \
  double p2k8, p2k9, p2k10, p2k11, p2k12, p2k13, p2k14;              \
                                                                     \
  double fac = A0[0];                                                \
  double *temp_ptr;                                                  \
  double *myA0, *myAnext;                                            \
                                                                     \
  int neg_x_slope, neg_y_slope, neg_z_slope;                         \
  int pos_x_slope, pos_y_slope, pos_z_slope;                         \
  int blockMin_x, blockMin_y, blockMin_z;                            \
  int blockMax_x, blockMax_y, blockMax_z;                            \
  int ii, jj, kk, i, j, k, t;                                        \
  int txx, tyy, tzz;                                                 \
                                                                     \
  for (kk= LENGTH; kk < nz-LENGTH; kk+= tz) {                        \
    neg_z_slope = LENGTH;                                            \
    pos_z_slope = -LENGTH;                                           \
    tzz = MIN(tz, nz-kk-LENGTH);                                     \
                                                                     \
    if (kk == LENGTH) {                                              \
      neg_z_slope = 0;                                               \
    }                                                                \
    if (kk == nz-tzz-LENGTH) {                                       \
      pos_z_slope = 0;                                               \
    }                                                                \
    for (jj= LENGTH; jj < ny-LENGTH; jj+= ty) {                      \
      neg_y_slope = LENGTH;                                          \
      pos_y_slope = -LENGTH;                                         \
      tyy = MIN(ty, ny-jj-LENGTH);                                   \
                                                                     \
      if (jj == LENGTH) {                                            \
        neg_y_slope = 0;                                             \
      }                                                              \
      if (jj == ny-tyy-LENGTH) {                                     \
        pos_y_slope = 0;                                             \
      }                                                              \
      for (ii= LENGTH + PAD(LENGTH); ii < nx-LENGTH - PAD(LENGTH); ii+= tx) { \
        neg_x_slope = LENGTH;                                        \
        pos_x_slope = -LENGTH;                                       \
        txx = MIN(tx, nx-ii-LENGTH - PAD(LENGTH));                   \
                                                                     \
        if (ii == LENGTH + PAD(LENGTH)) {                            \
          neg_x_slope = 0;                                           \
        }                                                            \
        if (ii == nx-txx-LENGTH - PAD(LENGTH)) {                     \
          pos_x_slope = 0;                                           \
        }                                                            \
                                                                     \
        myA0 = A0;                                                   \
        myAnext = Anext;                                             \
                                                                     \
        for (t= 0; t < timesteps; t++) {                             \
          blockMin_x = MAX(LENGTH + PAD(LENGTH), ii - t * neg_x_slope); \
          blockMin_y = MAX(LENGTH, jj - t * neg_y_slope);            \
          blockMin_z = MAX(LENGTH, kk - t * neg_z_slope);            \
                                                                     \
          blockMax_x = MAX(LENGTH + PAD(LENGTH), ii + txx + t * pos_x_slope); \
          blockMax_y = MAX(LENGTH, jj + tyy + t * pos_y_slope);      \
          blockMax_z = MAX(LENGTH, kk + tzz + t * pos_z_slope);      \
                                                                     \
          SEMISTENCIL_( myAnext, myA0, LENGTH,                       \
                        blockMin_x, blockMax_x,                      \
                        blockMin_y, blockMax_y,                      \
                        blockMin_z, blockMax_z )                     \
                                                                     \
          temp_ptr = myA0;                                           \
          myA0 = myAnext;                                            \
          myAnext = temp_ptr;                                        \
        }                                                            \
      }                                                              \
    }                                                                \
  }

#endif


/*
  StencilProbe Heat Equation
  Implements 7pt Timeskew Semi-stencil (1)
*/
void StencilProbe_timeskew_semi_1(double *A0, double *Anext, int nx, int ny, int nz,
                                  int tx, int ty, int tz, int timesteps) {

#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  TIMESKEW_SEMI( 1 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 13pt Timeskew Semi-stencil (2)
*/
void StencilProbe_timeskew_semi_2(double *A0, double *Anext, int nx, int ny, int nz,
                                  int tx, int ty, int tz, int timesteps) {

#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  TIMESKEW_SEMI( 2 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 25pt Timeskew Semi-stencil (4)
*/
void StencilProbe_timeskew_semi_4(double *A0, double *Anext, int nx, int ny, int nz,
                                  int tx, int ty, int tz, int timesteps) {

  TIMESKEW_SEMI( 4 );
}


/*
  StencilProbe Heat Equation
  Implements 43pt Timeskew Semi-stencil (7)
*/
void StencilProbe_timeskew_semi_7(double *A0, double *Anext, int nx, int ny, int nz,
                                  int tx, int ty, int tz, int timesteps) {

#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  TIMESKEW_SEMI( 7 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 85pt Timeskew Semi-stencil (14)
*/
void StencilProbe_timeskew_semi_14(double *A0, double *Anext, int nx, int ny, int nz,
                                   int tx, int ty, int tz, int timesteps) {

#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  TIMESKEW_SEMI( 14 );
#endif
}


/* This method traverses all of the cache blocks in a specific order to preserve
   dependencies.  For each cache block, it performs (possibly) several iterations while
   still respecting boundary conditions.
   NOTE: Positive slopes indicate that each iteration goes further out from the center
   of the current cache block, while negative slopes go toward the block center. */
#ifdef STENCILTEST
void StencilProbe_timeskew_semi(double *A0, double *Anext, int nx, int ny, int nz,
                           int tx, int ty, int tz, int timesteps, int length) {
#else
void StencilProbe(double *A0, double *Anext, int nx, int ny, int nz,
                  int tx, int ty, int tz, int timesteps, int length) {
#endif
  switch(length)
  {
    case 1: StencilProbe_timeskew_semi_1(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 2: StencilProbe_timeskew_semi_2(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 4: StencilProbe_timeskew_semi_4(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 7: StencilProbe_timeskew_semi_7(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    case 14:StencilProbe_timeskew_semi_14(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);break;
    default: printf("ABORTING: StencilProbe_timeskew_semi_%d NOT IMPLEMENTED!\n", length); exit(-1);
  }
}

