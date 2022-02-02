//
// standard header files
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//
// CUDA header file
//

#include <cuda.h>

//
// my header file
//

#include "poissinv_cuda.h"

//
// function prototype for quad precision evaluation of Poisson CDF
//

extern "C" void poissinv_quad(int, float, float*, float*, double*, double*);

//
// CUDA kernels
//

__global__ void poissinvf_bisection( int N, float lam,
                                     float *ulo_d, float *uhi_d ) {

  float x, xt, u_lo, u_hi, u_mid;
  int   tid = threadIdx.x + blockIdx.x*blockDim.x;

  if (tid < N) {
    u_hi  = 1.0f;
    u_lo  = 0.0f;
    u_mid = 0.5f*(u_hi + u_lo);
    xt    = (float) tid;

    while (u_mid>u_lo & u_mid<u_hi) {
      x = poissinvf(u_mid, lam);

      if (x>xt)
        u_hi = u_mid;
      else
        u_lo = u_mid;

      u_mid = 0.5f*(u_hi + u_lo);
    }
    ulo_d[tid] = u_lo;
    uhi_d[tid] = u_hi;
  }
}


__global__ void poissinv_bisection( int N, double lam,
                                    double *ulo_d, double *uhi_d ) {

  double x, xt, u_lo, u_hi, u_mid;
  int   tid = threadIdx.x + blockIdx.x*blockDim.x;

  if (tid < N) {
    u_hi  = 1.0;
    u_lo  = 0.0;
    u_mid = 0.5*(u_hi + u_lo);
    xt    = (double) tid;

    while (u_mid>u_lo & u_mid<u_hi) {
      x = poissinv(u_mid, lam);

      if (x>xt)
        u_hi = u_mid;
      else
        u_lo = u_mid;

      u_mid = 0.5*(u_hi + u_lo);
    }
    ulo_d[tid] = u_lo;
    uhi_d[tid] = u_hi;
  }
}


//////////////////////////////////////////////////
// main code
//////////////////////////////////////////////////

int main(int argc, char **argv) {

  float   lam;
  float  *ulo_h, *uhi_h, *ulo_ex, *uhi_ex;
  float  *ulo_d, *uhi_d;

  double *Ulo_h, *Uhi_h, *Ulo_ex, *Uhi_ex;
  double *Ulo_d, *Uhi_d;

  double err1;
  int    N, nblocks, nthreads; 

  // allocate memory

  int Nmax = 2000050;   // big enough for lambda up to 10^6

  ulo_ex = (float *)malloc(Nmax*sizeof(float));
  uhi_ex = (float *)malloc(Nmax*sizeof(float));
  Ulo_ex = (double *)malloc(Nmax*sizeof(double));
  Uhi_ex = (double *)malloc(Nmax*sizeof(double));

  ulo_h = (float *)malloc(Nmax*sizeof(float));
  uhi_h = (float *)malloc(Nmax*sizeof(float));
  cudaMalloc((void **)&ulo_d, Nmax*sizeof(float));
  cudaMalloc((void **)&uhi_d, Nmax*sizeof(float));

  Ulo_h = (double *)malloc(Nmax*sizeof(double));
  Uhi_h = (double *)malloc(Nmax*sizeof(double));
  cudaMalloc((void **)&Ulo_d, Nmax*sizeof(double));
  cudaMalloc((void **)&Uhi_d, Nmax*sizeof(double));

  // set values to test

  printf("       lam     SP_err      DP_err \n");

  lam = 0.5f;

  for (int count=0; count<20; count++) {
    lam = 2.0f*lam;
    N   = 50 + (int) (2*lam);
    if (N>Nmax) exit(1);

    // set number of blocks, and threads per block

    nthreads = 256;
    nblocks  = (N-1)/nthreads + 1;

//////////////////////////////////////////////////
// compute reference solution in quad precision
//////////////////////////////////////////////////

    poissinv_quad(N, lam, ulo_ex, uhi_ex, Ulo_ex, Uhi_ex);

//////////////////////////////////////////////////
// first do tests in single precision
//////////////////////////////////////////////////

    poissinvf_bisection<<<nblocks,nthreads>>>(N, lam, ulo_d, uhi_d);

    cudaMemcpy(ulo_h,ulo_d, N*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(uhi_h,uhi_d, N*sizeof(float),cudaMemcpyDeviceToHost);
    err1 = 0.0;

    for (int n=0; n<N; n++) {
      err1 += 0.5*fabs( (ulo_ex[n]-ulo_h[n])
                      + (uhi_ex[n]-uhi_h[n]));

      if (n>0) {
        if ( uhi_h[n] <= ulo_ex[n-1] ) {
          printf(" error: n = %d, uhi_h[n] = %20.16g, ulo_ex[n-1] = %20.16g \n",
                          n,      uhi_h[n],           ulo_ex[n-1]);
          exit(1);
        }
        if ( uhi_ex[n] <= ulo_h[n-1] ) {
          printf(" error: n = %d, uhi_ex[n] = %20.16g, ulo_h[n-1] = %20.16g \n",
                          n,      uhi_ex[n],           ulo_h[n-1]);
          exit(1);
        }
      }
    }

    printf("%10.4g   %9.3g   ",lam,err1);

//////////////////////////////////////////////////
// now re-do tests in double precision
//////////////////////////////////////////////////

    poissinv_bisection<<<nblocks,nthreads>>>(N, lam, Ulo_d, Uhi_d);

    cudaMemcpy(Ulo_h,Ulo_d, N*sizeof(double),cudaMemcpyDeviceToHost);
    cudaMemcpy(Uhi_h,Uhi_d, N*sizeof(double),cudaMemcpyDeviceToHost);

    err1 = 0.0;

    for (int n=0; n<N; n++) {
      err1 += 0.5*fabs( (Ulo_ex[n]-Ulo_h[n]) + (Uhi_ex[n]-Uhi_h[n]));

      if (n>0) {
        if ( Uhi_h[n] <= Ulo_ex[n-1] ) {
          printf(" error: n = %d, Uhi_h[n] = %20.16g, Ulo_ex[n-1] = %20.16g \n",
                          n,      Uhi_h[n],           Ulo_ex[n-1]);
          exit(1);
        }
        if ( Uhi_ex[n] <= Ulo_h[n-1] ) {
          printf(" error: n = %d, Uhi_ex[n] = %20.16g, Ulo_h[n-1] = %20.16g \n",
                          n,      Uhi_ex[n],           Ulo_h[n-1]);
          exit(1);
        }
      }
    }

    printf("%9.3g \n",err1);
  }

// free memory 

  cudaFree(Ulo_d);
  cudaFree(Uhi_d);
  free(Ulo_h);
  free(Uhi_h);

  cudaFree(ulo_d);
  cudaFree(uhi_d);
  free(ulo_h);
  free(uhi_h);

  free(ulo_ex);
  free(uhi_ex);
  free(Ulo_ex);
  free(Uhi_ex);

// CUDA exit -- needed to flush printf write buffer

  cudaDeviceReset();
}
