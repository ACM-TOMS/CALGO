/*
  Stencil Probe Testing version
  Version designed to test the resuls of each algorithm.
  Relative error (TOLERANCE) is used to check results.
*/

#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "util.h"
#include "cycle.h"
#ifdef HAVE_PAPI
#include <papi.h>
#endif
/* run.h has the run parameters */
#include "run.h"


int main(int argc,char *argv[])
{
  double *A0_naive, *A0_test;
  double *Anext_naive, *Anext_test;
  double *Afinal_naive, *Afinal_test;
  int nx, ny, nz;
  int tx, ty, tz;
  int timesteps, length;

  ticks t1, t2;
  double spt;
  
  /* parse command line options */
  if (argc < 9) {
    printf("\nUSAGE:\n%s <grid x> <grid y> <grid z> <block x> <block y> <block z> <timesteps> <length>\n", argv[0]);
    printf("\nLENGTH:\n 1-2-4-7-14.\n");
    printf("\nTIME SKEWING CONSTRAINTS:\nIn each dimension, <grid size - 2> should be a multiple of <block size>.\n");
    printf("\nCACHE OBLIVIOUS:\nblock x parameter is passed as CUTOFF in elements to the algorithm.\n");
    printf("\nCIRCULAR QUEUE CONSTRAINTS:\n<grid y - 2> should be a multiple of <block y>.  The block sizes in the other dimensions are ignored.\n\n");
    return EXIT_FAILURE;
  }

  length = atoi(argv[8]);

  /* add some padding for vectorial runs if needed (head/tail x) */
#ifdef VECSIZE
  if ((atoi(argv[1]) % VECSIZE) != 0) {
    printf("Error: X dimension (%d) is not multiple of VECSIZE (%d)\n", atoi(argv[1]), VECSIZE);
    return EXIT_FAILURE;
  }
#endif
  nx = atoi(argv[1]) + 2*length + 2*PAD(length);
  ny = atoi(argv[2]) + 2*length;
  nz = atoi(argv[3]) + 2*length;
  tx = atoi(argv[4]);
  ty = atoi(argv[5]);
  tz = atoi(argv[6]);
  timesteps = atoi(argv[7]);
  printf("[%s] %s: %dx%dx%d, blocking: %dx%dx%d, timesteps: %d, length: %d\n", VERSION, argv[0], nx, ny, nz, tx, ty, tz, timesteps, length);

#ifdef HAVE_PAPI
  PAPI_library_init(PAPI_VER_CURRENT);
# ifdef VECSIZE
  printf("*** padding: %d, VECSIZE: %d, nx: %d\n", PAD(length), VECSIZE, nx);
# endif
#endif
  
  /* find conversion factor from ticks to seconds */
  spt = seconds_per_tick();
  
  // allocate arrays
#ifdef __INTEL_COMPILER
  A0_naive=(double*)_mm_malloc(sizeof(double)*nx*ny*nz, 4096);
  A0_test=(double*)_mm_malloc(sizeof(double)*nx*ny*nz, 4096);
  Anext_naive=(double*)_mm_malloc(sizeof(double)*nx*ny*nz, 4096);
  Anext_test=(double*)_mm_malloc(sizeof(double)*nx*ny*nz, 4096);
#else
  A0_naive=(double*)malloc(sizeof(double)*nx*ny*nz);
  A0_test=(double*)malloc(sizeof(double)*nx*ny*nz);
  Anext_naive=(double*)malloc(sizeof(double)*nx*ny*nz);
  Anext_test=(double*)malloc(sizeof(double)*nx*ny*nz);
#endif

  // Run Naive Code
  StencilInit(nx,ny,nz,A0_naive);
  StencilInit(nx,ny,nz,Anext_naive);
  StencilProbe_naive(A0_naive, Anext_naive, nx, ny, nz, tx, ty, tz, timesteps, length);
  if (timesteps%2 == 0) {
    Afinal_naive = A0_naive;
  }
  else {
    Afinal_naive = Anext_naive;
  }
  
  printf("USING TIMER: %s \t  SECONDS PER TICK:%g \n", TIMER_DESC, spt);

  // Test Rivera Blocking
  StencilInit(nx,ny,nz,A0_test);
  StencilInit(nx,ny,nz,Anext_test);
  printf("Checking Rivera blocking...\n");
  StencilProbe_rivera(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps, length);
  if (timesteps%2 == 0) {
    Afinal_test = A0_test;
  }
  else {
    Afinal_test = Anext_test;
  }
  check_vals(Afinal_naive, Afinal_test, nx, ny, nz);

  // Test Semi-stencil
  StencilInit(nx,ny,nz,A0_test);
  StencilInit(nx,ny,nz,Anext_test);
  printf("Checking Semi-stencil...\n");
  StencilProbe_semi(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps, length);
  if (timesteps%2 == 0) {
    Afinal_test = A0_test;
  }
  else {
    Afinal_test = Anext_test;
  }
  check_vals(Afinal_naive, Afinal_test, nx, ny, nz);

  // Test Blocked Semi-stencil
  StencilInit(nx,ny,nz,A0_test);
  StencilInit(nx,ny,nz,Anext_test);
  printf("Checking Rivera Semi-stencil blocking...\n");
  StencilProbe_blocked_semi(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps, length);
  if (timesteps%2 == 0) {
    Afinal_test = A0_test;
  }
  else {
    Afinal_test = Anext_test;
  }
  check_vals(Afinal_naive, Afinal_test, nx, ny, nz);

  // Test Cache-Oblivious Semi-stencil Blocking
  StencilInit(nx,ny,nz,A0_test);
  StencilInit(nx,ny,nz,Anext_test);
  printf("Checking Cache-Oblivious Semi-stencil blocking...\n");
  StencilProbe_oblivious_semi(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps, length);
  if (timesteps%2 == 0) {
    Afinal_test = A0_test;
  }
  else {
    Afinal_test = Anext_test;
  }
  check_vals(Afinal_naive, Afinal_test, nx, ny, nz);

  // Test Time-Skewed Semi-stencil Blocking
  StencilInit(nx,ny,nz,A0_test);
  StencilInit(nx,ny,nz,Anext_test);
  printf("Checking Time-Skewed Semi-stencil blocking...\n");
  StencilProbe_timeskew_semi(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps, length);
  if (timesteps%2 == 0) {
    Afinal_test = A0_test;
  }
  else {
    Afinal_test = Anext_test;
  }
  check_vals(Afinal_naive, Afinal_test, nx, ny, nz);

  // Test Cache-Oblivious Blocking
  StencilInit(nx,ny,nz,A0_test);
  StencilInit(nx,ny,nz,Anext_test);
  printf("Checking Cache-Oblivious blocking...\n");
  StencilProbe_oblivious(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps, length);
  if (timesteps%2 == 0) {
    Afinal_test = A0_test;
  }
  else {
    Afinal_test = Anext_test;
  }
  check_vals(Afinal_naive, Afinal_test, nx, ny, nz);

  // Test Time-Skewed Blocking
  StencilInit(nx,ny,nz,A0_test);
  StencilInit(nx,ny,nz,Anext_test);
  printf("Checking Time-Skewed blocking...\n");
  StencilProbe_timeskew(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps, length);
  if (timesteps%2 == 0) {
    Afinal_test = A0_test;
  }
  else {
    Afinal_test = Anext_test;
  }
  check_vals(Afinal_naive, Afinal_test, nx, ny, nz);

// // Test Circular-Queue Blocking
// StencilInit(nx,ny,nz,A0_test);
// StencilInit(nx,ny,nz,Anext_test);
// if (timesteps > 1) {
//   CircularQueueInit(nx, ty, timesteps);
// }
// printf("Checking Circular-Queue blocking...\n");
// StencilProbe_circqueue(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps);
// Afinal_test = Anext_test;
// check_vals(Afinal_naive, Afinal_test, nx, ny, nz);
  
  /* free arrays */
#ifdef __INTEL_COMPILER
  _mm_free(Anext_naive);
  _mm_free(A0_naive);
  _mm_free(Anext_test);
  _mm_free(A0_test);
#else
  free(Anext_naive);
  free(A0_naive);
  free(Anext_test);
  free(A0_test);
#endif

  return EXIT_SUCCESS;
}

