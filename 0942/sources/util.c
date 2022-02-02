/*
  Stencil Probe utilities
  Helper functions for the probe.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "cycle.h"


#define MAX_SIZE 4*1024*1024  // 32 Mbytes (8 bytes elements)


/*
  This initializes the array A to be all 1's.  
  This is nearly superfluous (could use memset), but
  provides convenience and consistency nonetheless...
 */
void StencilInit(int nx, int ny, int nz, /* size of the array */
                 double *A) {            /* the array to initialize to 1s */
  long i;
  long last = nx*ny*nz;
  unsigned int seed = 1;

  /* Different parallel initializations should be considered depending on the
   * domain decomposition policy used among OpenMP threads IF MORE THAN ONE SOCKET IS USED!
   *
   * For instance:
   *
   * do i = 1, 4
   *  !$omp parallel num_threads(i)
   *    !$omp parallel num_threads(2)
   *      continue
   *    !$omp end parallel
   *  !$omp end parallel
   * end do
   *
   * #pragma omp parallel for private(i)
   */
#pragma omp parallel for default(shared) private(i, seed)
  for(i=0;i<last;i++) {
#ifdef RANDOMVALUES
    A[i]=(double)rand_r(&seed)/RAND_MAX;
#else
    A[i]=1.0;
#endif
  }
}

/*
  This function determines ticks per second.
  Inspired by OSKI function (bebop.cs.berkeley.edu)
*/
double seconds_per_tick()
{
  ticks t0,t1;
  unsigned int i = 3;
  double spt = 0;

  while (spt <= 0)
  {
    t0=getticks();  
    sleep(i);
    t1=getticks();
    spt = (double)i / elapsed(t1,t0);
    i++;
  }

  return spt;
}

/*
  Function to clear the cache, preventing data items in cache
  from making subsequent trials run faster.
*/
double clear_cache()
{
  int i;
  double* tarray, accum;

  /* Increase the buffer size if needed: now 32MB L3 */
  tarray = (double*) malloc(sizeof(double)*MAX_SIZE);

#pragma omp parallel for default(shared) private(i) reduction(+:accum)
  for (i=0; i<MAX_SIZE; i++) {
    tarray[i] += 1.0;
    accum += tarray[i];
  }

  free(tarray);

  /* We return accum to avoid compiler to change the tarray
   * initialization with a memset system call or be removed */
  return accum;
}

/*
  Handles error for PAPI routines
*/
void handle_error(int errcode)
{
  printf("Execution yielded error code %d : %s\n", errcode, ErrMsg[errcode]);
  exit(errcode);
}

/*
  Compare results between two executions using TOLERANCE as relative error
*/
void check_vals(double* A, double* B, int nx, int ny, int nz)
{
  int same, different;
  int i, j, k;
  double rel_err;

  same=different=0;

  for (k=0; k<nz; k++) {
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
        rel_err = fabs(A[Index3D(nx,ny,i,j,k)]-B[Index3D(nx,ny,i,j,k)])/A[Index3D(nx,ny,i,j,k)];
        if (rel_err < TOLERANCE)
          same++;
        else {
          different++;
          printf("at index %d %d %d --- A: %16e, B %16e\n", i, j, k,
                 A[Index3D(nx,ny,i,j,k)], B[Index3D(nx,ny,i,j,k)]);
        }
      }
    }
  }
  printf("Same: %d   Different: %d\n", same, different);
}

