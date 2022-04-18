/*
  Stencil Probe
  Main function.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "util.h"
#include "cycle.h"
#ifdef HAVE_PAPI
# include <papi.h>
#endif

/* run.h has the run parameters */
#include "run.h"


int main(int argc,char *argv[])
{
  int i;
  int nx, ny, nz;
  int tx, ty, tz;
  int timesteps, length;
  double *Anext;
  double *A0;
  
  ticks t1, t2;
  ticks tmin, tmax, tavg;
  ticks t[NUM_TRIALS];
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
#ifndef PLOT
  printf("[%s] %s: %dx%dx%d, blocking: %dx%dx%d, timesteps: %d, length: %d\n", VERSION, argv[0], nx, ny, nz, tx, ty, tz, timesteps, length);
# ifdef VECSIZE
  printf("*** padding: %d, VECSIZE: %d, nx: %d\n", PAD(length), VECSIZE, nx);
# endif
#endif

#ifdef TIMESKEW_BLK
  if ((tx < timesteps*length) || (ty < timesteps*length) || (tz < timesteps*length))
  {
    printf("Error in parameters for timeskewing: tx & ty & tz must be bigger or equal to timesteps*length\n");
    return EXIT_FAILURE;
  }
#endif

#ifdef HAVE_PAPI
  PAPI_library_init(PAPI_VER_CURRENT);
#endif
  
  /* find conversion factor from ticks to seconds */
  spt = seconds_per_tick();
  
  /* allocate arrays */ 
#ifdef __INTEL_COMPILER
  Anext = (double*)_mm_malloc(sizeof(double)*nx*ny*nz, 4096);
  A0    = (double*)_mm_malloc(sizeof(double)*nx*ny*nz, 4096);
#else
  Anext = (double*)malloc(sizeof(double)*nx*ny*nz);
  A0    = (double*)malloc(sizeof(double)*nx*ny*nz);
#endif
  
#ifndef PLOT
  printf("USING TIMER: %s \t  SECONDS PER TICK:%g \n", TIMER_DESC, spt);
#endif
  
  for (i=0;i<NUM_TRIALS;i++) {

#ifndef PEAK
    /* Initialize arrays to all ones
     * This initialization can be dangerous in
     * multi-socket executions, check OpenMP policy */
    StencilInit(nx,ny,nz,Anext);
    StencilInit(nx,ny,nz,A0);
# ifdef CIRCULARQUEUEPROBE
    //if (timesteps > 1) {                                                                                                                      
    //  CircularQueueInit(nx, ty, timesteps);                                                                                                   
    //}
# endif

    clear_cache();
#endif

    t1 = getticks();
 
    /* stencil function */ 
    StencilProbe(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps, length);
 
    t2 = getticks();
 
#ifndef PLOT
    printf("elapsed ticks: %g  time:%g \n", elapsed(t2, t1), spt * elapsed(t2,t1));
#else
    t[i] = elapsed(t2, t1);
#endif
  }

#ifdef PLOT
  /* Find min, max and average */
  tmin = t[0]; /* To avoid problems with different     */
  tmax = t[0]; /* types of ticks datatype in tmin/tmax */
  tavg = t[0];
  for (i=1;i<NUM_TRIALS;i++)
  {
    if (t[i] < tmin) tmin = t[i];
    if (t[i] > tmax) tmax = t[i];
    tavg += t[i];
  }

  /* BinName  VolSize (x,y,z)  BlkSize (tx,ty,tz) Ts Lg Mi Ma Av */
  printf( "%s\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d %14.8f %14.8f %14.8f\n",
          argv[0], nx-2*length-2*PAD(length), ny-2*length, nz-2*length, tx, ty, tz,
          timesteps, length, spt * tmin, spt * tmax, (spt * tavg)/NUM_TRIALS );
#endif

  /* free arrays */
#ifdef __INTEL_COMPILER
  _mm_free(Anext);
  _mm_free(A0);
#else
  free(Anext);
  free(A0);
#endif

  return EXIT_SUCCESS;
}

