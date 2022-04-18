/*
  Stencil Probe HWC version
  Version designed to run with HWC monitoring though PAPI, PAPIEX or EXTRAE in BG/P.
    PAPI tracing: -DPAPITRACE
    PAPIEX tracing: NOTHING
    EXTRAE tracing: -DSEQTRACE
  Set PAPI_COUNTERS environment variables with PAPI presets in a comma-separated fashion.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
  double *Anext;
  double *A0;
  int nx, ny, nz;
  int tx, ty, tz;
  int timesteps, length;

  TRACE_INIT()

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

  /* allocate arrays */
#ifdef __INTEL_COMPILER
  Anext=(double*)_mm_malloc(sizeof(double)*nx*ny*nz, 4096);
  A0=(double*)_mm_malloc(sizeof(double)*nx*ny*nz, 4096);
#else
  Anext=(double*)malloc(sizeof(double)*nx*ny*nz);
  A0=(double*)malloc(sizeof(double)*nx*ny*nz);
#endif

  TRACE_START();

  /* stencil function */ 
  StencilProbe(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps, length);

  TRACE_STOP();

  /* free arrays */
#ifdef __INTEL_COMPILER
  _mm_free(Anext);
  _mm_free(A0);
#else
  free(Anext);
  free(A0);
#endif

  TRACE_FINI()

  return EXIT_SUCCESS;
}

