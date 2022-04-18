//
// standard header files
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//
// linux timing routine
//

#include <sys/time.h>

inline double elapsed_time(double *et) {
  struct timeval t;

  double old_time = *et;

  gettimeofday( &t, (struct timezone *)0 );
  *et = t.tv_sec + t.tv_usec*1.0e-6;

  return *et - old_time;
}

//
// my header file
//

#include "poissinv.h"

void poissinv_test_scalar(int N, double lam) {
  double x, u;
  int    n;

  for (n=0; n<N; n++) {
    u = (n+0.5) / N;

    if (lam>400.0)
      x = normcdfinv_as241(u);
    else
      // extra n*1e-100 is to prevent compiler
      // optimising for fixed lam
      x = poissinv(u, lam+n*1e-100);

// needed to prevent compiler discarding everything
    if (x==-999.0) printf("negative x\n");
  }
}

void poissinv_test_vector(int N, double lam) {
  double x, u;
  int    n;

  for (n=0; n<N; n++) {
    u = (n+0.5) / N;

    if (lam>400.0)
      x = normcdfinv_as241(u);
    else
      // extra n*1e-100 is to prevent compiler
      // optimising for fixed lam
      x = poissinv_v(u, lam+n*1e-100);

// needed to prevent compiler discarding everything
    if (x==-999.0) printf("negative x\n");
  }
}



//
// main code
//

int main(int argc, char **argv) {
  double timer, elapsed;  // timer variable and elapsed time
  double lam;
  int    N, pass, count, Count=6; 

  N = (1<<24);

  // execute code

  for (pass=0; pass<2; pass++) {
    lam = 0.125;

    if (pass==0)
      printf("\nscalar algorithm performance tests \n");
    else
      printf("\nvector algorithm performance tests \n");
    printf("---------------------------------- \n");
    printf("  lambda   execution time   samples/sec \n");

    for (count=0; count<Count; count++) {
      lam = lam*4.0;

      elapsed_time(&timer);  // initialise timer

      if (pass==0)
        poissinv_test_scalar(N, lam);
      else
        poissinv_test_vector(N, lam);

      elapsed = elapsed_time(&timer);
 
      if (count==Count-1)
        printf(" normcdfinv  %9.4f     %10.3g \n",
                   elapsed, N/elapsed);
      else if (count>0)  // skip first one
        printf("   %4g      %9.4f     %10.3g \n",
              lam, elapsed, N/elapsed);
    }
  }

  return 0;
}
