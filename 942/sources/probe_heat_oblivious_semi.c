/*
  StencilProbe Heat Equation (Cache Oblivious + Semi-stencil version)
  Implements stencil from Chombo's heattut example with oblivious & semi-stencil.
  Raúl de la Cruz (delacruz@bsc.es)
  Barcelona Supercomputing Center

  Based on Cache Oblivious algorithm from "Cache Oblivious Stencil Computations"
  of Matteo Frigo and Volker Strumpen + CUTOFF (stencilprobe feature).
  This version never cuts over unit-stride dimension X.
  OpenMP version: based on "The Cache Complexity of Multithreaded Cache Oblivious Algorithms"
  paper of Matteo Frigo and Volker Strumpen. Parallelism is created in space cuts (either
  Z or Y dimensions depending on the bigger dimension), and two types of trapezoids
  (black and grey) are computed by each thread to preserve dependencies.
  Nested parallelism is not allowed (parallelism is open only once).

  Trapezoid space cut: 2 balanced trapezoids (T1 and T2) are created (same number of points
  to compute) through the center point (c(t,x), depicted as x in cut examples of below).

    Trapezoid: T(t0,t1,x0,dx0,x1,dx1) where t0 <= t1, x0 <= x1 and x0+dx0*dt <= x1+dx1*dt
    Height: dt (t1 - t0)
    Width: (x1-x0) + (dx1-dx0)*dt/2
    Center: c(t,x) where t = dt/2 and x = (x0+x1)/2 + (dx0+dx1)*dt/4
    zm (forward projection of center point): (2*(z0+z1) + (2*ds+dz0+dz1)*dt)/4
    zm' (backward projection of center point): (2*(z0+z1) + (2*-ds+dz0+dz1)*dt)/4

  To ensure the creation of two well-defined trapezoids, the four edges of each new trapezoid
  must fullfil the following 4 conditions:

    z0 < zm, zm < z1, z0+dx0*dt < zm' and zm < z1+dx1*dt
    where T1(t0,t1,z0,dx0,zm,-ds) and T2(t0,t1,zm,-ds,z1,dx1)

    * Examples of well-defined trapezoid cuts

                zm'                                   zm'
      z0+dx0*dt ·-\----· z1+dx1*dt      z0+dx0*dt ·-----\---------· z1+dx1*dt
               /   \    \                          \     \  T2   /
              / T1  x T2 \                          \     x     /
             /       \    \                          \  T1 \   /
            ·---------\----·                          \_____\_/
           z0          zm   z1                       z0    zm\ z1

    * Examples of wrong trapezoid cuts

      zm'  z0+dx0*dt                           zm'
          \·-· z1+dx1*dt         z0+dx0*dt ·----\--------· z1+dx1*dt
          /\T2\                             \    \  T2  /
         /  x  \                             \    x    /
        / T1 \  \                             \ T1 \  /
       ·------\--·                             \    \/
      z0      zm  z1                            \___/\ zm
                                               z0   z1
*/
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "semi.h"


#ifdef _OPENMP

#define OBLIVIOUS_SEMI( WALK, CUTOFF, ds )                                                  \
  double p2c;                                                                               \
  double p2i1, p2i2, p2i3,  p2i4,  p2i5,  p2i6,  p2i7;                                      \
  double p2i8, p2i9, p2i10, p2i11, p2i12, p2i13, p2i14;                                     \
  double p2j1, p2j2, p2j3,  p2j4,  p2j5,  p2j6,  p2j7;                                      \
  double p2j8, p2j9, p2j10, p2j11, p2j12, p2j13, p2j14;                                     \
  double p2k1, p2k2, p2k3,  p2k4,  p2k5,  p2k6,  p2k7;                                      \
  double p2k8, p2k9, p2k10, p2k11, p2k12, p2k13, p2k14;                                     \
                                                                                            \
  int dt = t1-t0;                                                                           \
  int wx = ((x1-x0) + (dx1-dx0) * dt * 0.5);          /* Compute 3D trapezoid volume */     \
  int wy = ((y1-y0) + (dy1-dy0) * dt * 0.5);          /* for CUTOFF parameter and    */     \
  int wz = ((z1-z0) + (dz1-dz0) * dt * 0.5);          /* serial space cuts           */     \
  int vol = wx * wy * wz;                                                                   \
                                                                                            \
  if (dt > 1) {                                                                             \
    int r, l, rl, yy, zz;                                                                   \
    if (!omp_in_parallel()) r = omp_get_max_threads();  /* Parallelism is open only once */ \
    else r = 1;                                                                             \
                                                                                            \
    if (r != 1) {                                                /* Parallel code */        \
      if ( (z1-z0) >= 2 * ds * dt * r) {                         /* Parallel Space Z-cut */ \
                                                                                            \
        _Pragma("omp parallel default(shared) private(l,rl,zz)")                            \
        {                                                                                   \
          l = (z1-z0) / r;                                                                  \
          rl = (z1-z0) % r;                                                                 \
          if (rl > omp_get_thread_num()) {                                                  \
            l++;                                                                            \
            zz = z0 + l*omp_get_thread_num();                                               \
          }                                                                                 \
          else zz = z0 + l*omp_get_thread_num() + rl;                                       \
                                                                                            \
          WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,          /* Spawn black trapezoids */ \
               y0,dy0,y1,dy1,zz,ds,MIN(zz+l,z1),-ds);                                       \
                                                         /* Wait black trap. to complete */ \
          _Pragma("omp barrier")                         /* before spawning grey trapez. */ \
                                                                                            \
                                                               /* Spawn grey trapezoids  */ \
          if (omp_get_thread_num() == 0) {                     /* (2 halves trap. th. 0) */ \
            WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,y0,dy0,y1,dy1,z0,dz0,z0,ds);         \
            WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,y0,dy0,y1,dy1,z1,-ds,z1,dz1);        \
          }                                                                                 \
          else {                                                                            \
            WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,y0,dy0,y1,dy1,zz,-ds,zz,ds);         \
          }                                                                                 \
        }                                                                                   \
      }                                                                                     \
      else if ( (y1-y0) >= 2 * ds * dt * r) {                    /* Parallel Space Y-cut */ \
                                                                                            \
        _Pragma("omp parallel default(shared) private(l,rl,yy)")                            \
        {                                                                                   \
          l = (y1-y0) / r;                                                                  \
          rl = (y1-y0) % r;                                                                 \
          if (rl > omp_get_thread_num()) {                                                  \
            l++;                                                                            \
            yy = y0 + l*omp_get_thread_num();                                               \
          }                                                                                 \
          else yy = y0 + l*omp_get_thread_num() + rl;                                       \
                                                                                            \
          WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,          /* Spawn black trapezoids */ \
               yy,ds,MIN(yy+l,y1),-ds,z0,dz0,z1,dz1);                                       \
                                                         /* Wait black trap. to complete */ \
          _Pragma("omp barrier")                         /* before spawning grey trapez. */ \
                                                                                            \
                                                               /* Spawn grey trapezoids  */ \
          if (omp_get_thread_num() == 0) {                     /* (2 halves trap. th. 0) */ \
            WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,y0,dy0,y0,ds,z0,dz0,z1,dz1);         \
            WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,y1,-ds,y1,dy1,z0,dz0,z1,dz1);        \
          }                                                                                 \
          else {                                                                            \
            WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,yy,-ds,yy,ds,z0,dz0,z1,dz1);         \
          }                                                                                 \
        }                                                                                   \
      }                                                                                     \
      else {                                                     /* Time cut */             \
        int s = dt/2;                                                                       \
        WALK(A,nx,ny,nz,CUTOFF,t0,t0+s,x0,dx0,x1,dx1,y0,dy0,y1,dy1,z0,dz0,z1,dz1);          \
        WALK(A,nx,ny,nz,CUTOFF,t0+s,t1,x0+dx0*s,dx0,x1+dx1*s,dx1,y0+dy0*s,dy0,              \
             y1+dy1*s,dy1,z0+dz0*s,dz0,z1+dz1*s,dz1);                                       \
      }                                                                                     \
    }                                    /* Serial code, check well-defined. Cut black */   \
    else {                               /* (longer base) & grey (smaller base) trapez.*/   \
      if (wz >= 2 * ds * dt) {                                   /* Serial Space Z-cut */   \
        int zm  = (2 * (z0+z1) + (2*ds+dz0+dz1) * dt) / 4;                                  \
        WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,y0,dy0,y1,dy1,z0,dz0,zm,-ds);            \
        WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,y0,dy0,y1,dy1,zm,-ds,z1,dz1);            \
      }                                                                                     \
      else if (wy >= 2 * ds * dt) {                              /* Serial Space Y-cut */   \
        int ym  = (2 * (y0+y1) + (2*ds+dy0+dy1) * dt) / 4;                                  \
        WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,y0,dy0,ym,-ds,z0,dz0,z1,dz1);            \
        WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,ym,-ds,y1,dy1,z0,dz0,z1,dz1);            \
      }                                                                                     \
      else {                                                     /* Time cut */             \
        int s = dt/2;                                                                       \
        WALK(A,nx,ny,nz,CUTOFF,t0,t0+s,x0,dx0,x1,dx1,y0,dy0,y1,dy1,z0,dz0,z1,dz1);          \
        WALK(A,nx,ny,nz,CUTOFF,t0+s,t1,x0+dx0*s,dx0,x1+dx1*s,dx1,y0+dy0*s,dy0,              \
             y1+dy1*s,dy1,z0+dz0*s,dz0,z1+dz1*s,dz1);                                       \
      }                                                                                     \
    }                                                                                       \
  }                                                                                         \
  else if (dt == 1 || vol < CUTOFF) {                            /* Base case */            \
    int i, j, k, t;                                                                         \
    int xs, xe, ys, ye, zs, ze;                                                             \
    double fac = A[0][0];                                                                   \
    double *Anext, *A0;                                                                     \
    for (t= t0; t< t1; t++) {                                                               \
      Anext = A[(t+1)%2];                             /* Workarround for Intel compilers */ \
      A0    = A[t%2];                                 /* and weird behavior with perfor- */ \
      xs = x0+(t-t0)*dx0; xe = x1+(t-t0)*dx1;         /* mance on parameters in macros.  */ \
      ys = y0+(t-t0)*dy0; ye = y1+(t-t0)*dy1;                                               \
      zs = z0+(t-t0)*dz0; ze = z1+(t-t0)*dz1;                                               \
      SEMISTENCIL_( Anext, A0, ds,                                                          \
                    xs, xe, ys, ye, zs, ze )                                                \
    }                                                                                       \
  }

#else

#define OBLIVIOUS_SEMI( WALK, CUTOFF, ds )                                       \
  double p2c;                                                                    \
  double p2i1, p2i2, p2i3,  p2i4,  p2i5,  p2i6,  p2i7;                           \
  double p2i8, p2i9, p2i10, p2i11, p2i12, p2i13, p2i14;                          \
  double p2j1, p2j2, p2j3,  p2j4,  p2j5,  p2j6,  p2j7;                           \
  double p2j8, p2j9, p2j10, p2j11, p2j12, p2j13, p2j14;                          \
  double p2k1, p2k2, p2k3,  p2k4,  p2k5,  p2k6,  p2k7;                           \
  double p2k8, p2k9, p2k10, p2k11, p2k12, p2k13, p2k14;                          \
                                                                                 \
  int dt = t1-t0;                                                                \
  int wx = ((x1-x0) + (dx1-dx0) * dt * 0.5);  /* Compute 3D trapezoid volume */  \
  int wy = ((y1-y0) + (dy1-dy0) * dt * 0.5);  /* for CUTOFF parameter        */  \
  int wz = ((z1-z0) + (dz1-dz0) * dt * 0.5);                                     \
  int vol = wx * wy * wz;                                                        \
                                                                                 \
  if (dt == 1 || vol < CUTOFF) {                              /* Base case */    \
    int i, j, k, t;                                                              \
    int xs, xe, ys, ye, zs, ze;                                                  \
    double fac = A[0][0];                                                        \
    double *Anext, *A0;                                                          \
    for (t = t0; t<t1;t++) {                                                     \
      Anext = A[(t+1)%2];                    /* Workarround for Intel compilers*/\
      A0    = A[t%2];                        /* and weird behavior with perfor-*/\
      xs = x0+(t-t0)*dx0; xe = x1+(t-t0)*dx1;/* mance on parameters in macros. */\
      ys = y0+(t-t0)*dy0; ye = y1+(t-t0)*dy1;                                    \
      zs = z0+(t-t0)*dz0; ze = z1+(t-t0)*dz1;                                    \
      SEMISTENCIL_( Anext, A0, ds,                                               \
                    xs, xe, ys, ye, zs, ze )                                     \
    }                                                                            \
  }                                                                              \
  else if (dt > 1) {                                                             \
    if (wz >= 2 * ds * dt) {                                  /* Space Z-cut */  \
      int zm = (2* (z0+z1) + (2*ds+dz0+dz1) * dt) / 4;                           \
      WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,y0,dy0,y1,dy1,z0,dz0,zm,-ds);   \
      WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,y0,dy0,y1,dy1,zm,-ds,z1,dz1);   \
    }                                                                            \
    else if (wy >= 2 * ds * dt) {                             /* Space Y-cut */  \
      int ym = (2* (y0+y1) + (2*ds+dy0+dy1) * dt) / 4;                           \
      WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,y0,dy0,ym,-ds,z0,dz0,z1,dz1);   \
      WALK(A,nx,ny,nz,CUTOFF,t0,t1,x0,dx0,x1,dx1,ym,-ds,y1,dy1,z0,dz0,z1,dz1);   \
    }                                                                            \
    else {                                                    /* Time cut */     \
      int s = dt/2;                                                              \
      WALK(A,nx,ny,nz,CUTOFF,t0,t0+s,x0,dx0,x1,dx1,y0,dy0,y1,dy1,z0,dz0,z1,dz1); \
      WALK(A,nx,ny,nz,CUTOFF,t0+s,t1,x0+dx0*s,dx0,x1+dx1*s,dx1,y0+dy0*s,dy0,     \
           y1+dy1*s,dy1,z0+dz0*s,dz0,z1+dz1*s,dz1);                              \
    }                                                                            \
  }

#endif


/*
  StencilProbe Heat Equation
  Implements 7pt Oblivious Semi-stencil (1)
*/
void semi_walk3_1(double* A[], int nx, int ny, int nz, int CUTOFF,
                  int t0, int t1, int x0, int dx0, int x1, int dx1,
                  int y0, int dy0, int y1, int dy1,
                  int z0, int dz0, int z1, int dz1) {

#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  OBLIVIOUS_SEMI( semi_walk3_1, CUTOFF, 1 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 13pt Oblivious Semi-stencil (2)
*/
void semi_walk3_2(double* A[], int nx, int ny, int nz, int CUTOFF,
                  int t0, int t1, int x0, int dx0, int x1, int dx1,
                  int y0, int dy0, int y1, int dy1,
                  int z0, int dz0, int z1, int dz1) {

#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  OBLIVIOUS_SEMI( semi_walk3_2, CUTOFF, 2 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 25pt Oblivious Semi-stencil (4)
*/
void semi_walk3_4(double* A[], int nx, int ny, int nz, int CUTOFF,
                  int t0, int t1, int x0, int dx0, int x1, int dx1,
                  int y0, int dy0, int y1, int dy1,
                  int z0, int dz0, int z1, int dz1) {

  OBLIVIOUS_SEMI( semi_walk3_4, CUTOFF, 4 );
}


/*
  StencilProbe Heat Equation
  Implements 43pt Oblivious Semi-stencil (7)
*/
void semi_walk3_7(double* A[], int nx, int ny, int nz, int CUTOFF,
                  int t0, int t1, int x0, int dx0, int x1, int dx1,
                  int y0, int dy0, int y1, int dy1,
                  int z0, int dz0, int z1, int dz1) {

#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  OBLIVIOUS_SEMI( semi_walk3_7, CUTOFF, 7 );
#endif
}


/*
  StencilProbe Heat Equation
  Implements 85pt Oblivious Semi-stencil (14)
*/
void semi_walk3_14(double* A[], int nx, int ny, int nz, int CUTOFF,
                   int t0, int t1, int x0, int dx0, int x1, int dx1,
                   int y0, int dy0, int y1, int dy1,
                   int z0, int dz0, int z1, int dz1) {

#if (!defined(SSE) && !defined(AVX) && !defined(MIC))
  OBLIVIOUS_SEMI( semi_walk3_14, CUTOFF, 14 );
#endif
}


/*
 * NOTE: tx parameter is used as CUTOFF in Cache Oblivious algorithm
 */

#ifdef STENCILTEST
void StencilProbe_oblivious_semi(double* A0, double* Anext, int nx, int ny, int nz,
                                 int tx, int ty, int tz, int timesteps, int length) {
#else
void StencilProbe(double *A0, double *Anext, int nx, int ny, int nz,
                  int tx, int ty, int tz, int timesteps, int length) {
#endif
  double* A[2] = {A0, Anext};

  switch(length)
  {
    case 1: semi_walk3_1(A, nx, ny, nz, tx,
                         0, timesteps,
                         1 + PAD(1), 0, nx-1 - PAD(1), 0,
                         1, 0, ny-1, 0,
                         1, 0, nz-1, 0);break;
    case 2: semi_walk3_2(A, nx, ny, nz, tx,
                         0, timesteps,
                         2 + PAD(2), 0, nx-2 - PAD(2), 0,
                         2, 0, ny-2, 0,
                         2, 0, nz-2, 0);break;
    case 4: semi_walk3_4(A, nx, ny, nz, tx,
                         0, timesteps,
                         4 + PAD(4), 0, nx-4 - PAD(4), 0,
                         4, 0, ny-4, 0,
                         4, 0, nz-4, 0);break;
    case 7: semi_walk3_7(A, nx, ny, nz, tx,
                         0, timesteps,
                         7 + PAD(7), 0, nx-7 - PAD(7), 0,
                         7, 0, ny-7, 0,
                         7, 0, nz-7, 0);break;
    case 14: semi_walk3_14(A, nx, ny, nz, tx,
                           0, timesteps,
                           14 + PAD(14), 0, nx-14 - PAD(14), 0,
                           14, 0, ny-14, 0,
                           14, 0, nz-14, 0);break;
    default: printf("ABORTING: StencilProbe_oblivious_semi_%d NOT IMPLEMENTED!\n", length); exit(-1);
  }
}

