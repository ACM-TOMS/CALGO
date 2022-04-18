//
// standard header files
//

#include <stdio.h>
#include <stdlib.h>

//
// CUDA header file
//

#include <cuda.h>

//
// my header file
//

#include "poissinv_cuda.h"

__global__ void poissinvf_test(int N, float lam) {

  float x, u;
  int   tid = threadIdx.x + blockIdx.x*blockDim.x;

  u = (tid + 0.5f) / N;

  for (int n=0; n<10; n++) {
    u += 1e-30f;

    if (lam>800.0f)
      x = normcdfinvf(u);
    else {
      if (lam>400.0f) {
        int n = 1 << (tid & 7); 
        lam   = 1.0f * (float) n;
      }
      x = poissinvf(u, lam);
    }

// needed to prevent compiler discarding everything
    if (x==-999.0f) printf("negative x\n");
  }
}

__global__ void poissinv_test(int N, float lam) {

  float x, u;
  int   tid = threadIdx.x + blockIdx.x*blockDim.x;

  u = (tid + 0.5f) / N;

  for (int n=0; n<10; n++) {
    u += 1e-30f;

    if (lam>800.0f)
      x = normcdfinv((double) u);
    else {
      if (lam>400.0f) {
        int n = 1 << (tid & 7); 
        lam   = 1.0f * (float) n;
      }
      x = poissinv((double) u, (double) lam);
    }

// needed to prevent compiler discarding everything
    if (x==-999.0f) printf("negative x\n");
  }
}



//
// main code
//

int main(int argc, char **argv) {
  float lam;
  int   N, nblocks, nthreads, Count=6; 

// CUDA timing

  float milli;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // set number of blocks, and threads per block

  N = (1<<24);
  nthreads = 256;
  nblocks  = N / nthreads;

  // execute kernels

  for (int pass=0; pass<2; pass++) {
    if (pass==0)
      printf("\nsingle precision performance tests \n");
    else
      printf("\ndouble precision performance tests \n");
    printf("---------------------------------- \n");
    printf("  lambda   execution time   samples/sec \n");

    lam = 0.125f;

    for (int count=0; count<=Count; count++) {
      lam = lam*4.0f;

      cudaEventRecord(start);

      if (pass==0)
        poissinvf_test<<<nblocks,nthreads>>>(N, lam);
      else
        poissinv_test<<<nblocks,nthreads>>>(N, lam);

      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milli, start, stop);

// factor 10 due to repeat in test routines
// factor 1e3 due to timing in milliseconds

      if (lam>0.5f) {        // skip first one for more accurate timing
        if (count==Count)
          printf("\n normcdfinv  %9.4f     %10.3g \n",
                     milli, float(N)*10.0*1e3/milli);
        else if (count==Count-1)
          printf("   mixed     %9.4f     %10.3g \n",
                     milli, float(N)*10.0*1e3/milli);
        else
          printf("   %4g      %9.4f     %10.3g \n",
                lam, milli, float(N)*10.0*1e3/milli);
      }
    }

  }

// CUDA exit -- needed to flush printf write buffer

  cudaDeviceReset();
  return 0;
}
