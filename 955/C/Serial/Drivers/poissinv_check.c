//
// standard header files
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//
// my header file
//

#include "poissinv.h"

//
// function prototype for quad precision evaluation of Poisson CDF
//

void poissinv_quad(int, float, float*, float*, double*, double*);
void poisscinv_quad(int, float, double*, double*);

//


void poissinv_bisection_scalar( int tid, int N, double lam,
                                    double *ulo_d, double *uhi_d ) {

  double x, xt, u_lo, u_hi, u_mid;

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

void poissinv_bisection_vector( int tid, int N, double lam,
                                    double *ulo_d, double *uhi_d ) {

  double x, xt, u_lo, u_hi, u_mid;

  if (tid < N) {
    u_hi  = 1.0;
    u_lo  = 0.0;
    u_mid = 0.5*(u_hi + u_lo);
    xt    = (double) tid;

    while (u_mid>u_lo & u_mid<u_hi) {
      x = poissinv_v(u_mid, lam);

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


void poisscinv_bisection_scalar( int tid, int N, double lam,
                                    double *ulo_d, double *uhi_d ) {

  double x, xt, u_lo, u_hi, u_mid;

  if (tid < N) {
    u_hi  = 1.0;
    u_lo  = 0.0;
    u_mid = 0.5*(u_hi + u_lo);
    xt    = (double) tid;

    while (u_mid>u_lo & u_mid<u_hi) {
      x = poisscinv(u_mid, lam);

      if (x<xt)
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
  float                  *ulo_ex, *uhi_ex;
  double  Lam;
  double *Ulo_h, *Uhi_h, *Ulo_ex, *Uhi_ex;

  double err1;
  int    N, n, count;

  // allocate memory

  int Nmax = 2001000;   // big enough for lambda up to 10^6

  ulo_ex = (float  *)malloc(Nmax*sizeof(float));
  uhi_ex = (float  *)malloc(Nmax*sizeof(float));
  Ulo_ex = (double *)malloc(Nmax*sizeof(double));
  Uhi_ex = (double *)malloc(Nmax*sizeof(double));
  Ulo_h  = (double *)malloc(Nmax*sizeof(double));
  Uhi_h  = (double *)malloc(Nmax*sizeof(double));

  // set values to test

  printf("       lam     scalar_err   vector_err \n");

  lam = 0.5f;

  for (count=0; count<20; count++) {
    lam = 2.0f*lam;
    N   = 50 + (int) (2*lam);
    if (N>Nmax) exit(1);

//////////////////////////////////////////////////
// compute reference solution in quad precision
//////////////////////////////////////////////////

    poissinv_quad(N, lam, ulo_ex, uhi_ex, Ulo_ex, Uhi_ex);

//////////////////////////////////////////////////
// first check scalar version
//////////////////////////////////////////////////

    err1 = 0.0;
    Lam  = lam;

    for (n=0; n<N; n++) {
      poissinv_bisection_scalar(n, N, Lam, Ulo_h, Uhi_h);

      err1 += 0.5*fabs( (Ulo_ex[n]-Ulo_h[n]) + (Uhi_ex[n]-Uhi_h[n]));

      if (n>0) {
        if ( Uhi_h[n] <= Ulo_ex[n-1] ) {
          printf("\n error: lam = %f, n = %d, Uhi_h[n] = %20.16g, Ulo_ex[n-1] = %20.16g \n",
                            lam,      n,      Uhi_h[n],           Ulo_ex[n-1]);
          exit(1);
        }
        if ( Uhi_ex[n] <= Ulo_h[n-1] ) {
          printf("\n error: lam = %f, n = %d, Uhi_ex[n] = %20.16g, Ulo_h[n-1] = %20.16g \n",
		            lam,      n,      Uhi_ex[n],           Ulo_h[n-1]);
          exit(1);
        }
      }
    }

    printf("%10.4g     %9.3g ",lam,err1);

//////////////////////////////////////////////////
// then check vector version
//////////////////////////////////////////////////

    err1 = 0.0;
    Lam  = lam;

    for (n=0; n<N; n++) {
      poissinv_bisection_vector(n, N, Lam, Ulo_h, Uhi_h);

      err1 += 0.5*fabs( (Ulo_ex[n]-Ulo_h[n]) + (Uhi_ex[n]-Uhi_h[n]));

      if (n>0) {
        if ( Uhi_h[n] <= Ulo_ex[n-1] ) {
          printf("\n error: n = %d, Uhi_h[n] = %20.16g, Ulo_ex[n-1] = %20.16g \n",
                            n,      Uhi_h[n],           Ulo_ex[n-1]);
          exit(1);
        }
        if ( Uhi_ex[n] <= Ulo_h[n-1] ) {
          printf("\n error: n = %d, Uhi_ex[n] = %20.16g, Ulo_h[n-1] = %20.16g \n",
                            n,      Uhi_ex[n],           Ulo_h[n-1]);
          exit(1);
        }
      }
    }

    printf("   %9.3g \n",err1);
  }


//////////////////////////////////////////////////
// re-do for complementary version
//////////////////////////////////////////////////

  // set values to test

  printf("\n       lam     complement_err \n");

  lam = 0.5f;

  for (count=0; count<20; count++) {
    lam = 2.0f*lam;
    N   = 1000 + (int) (2*lam);
    if (N>Nmax) exit(1);

//////////////////////////////////////////////////
// compute reference solution in quad precision
//////////////////////////////////////////////////

    poisscinv_quad(N, lam, Ulo_ex, Uhi_ex);

/*
    if (count==0)
      for (n=0; n<N; n++) {
        printf(" %g  %g  \n",Ulo_ex[n], Uhi_ex[n]);
      }
*/

//////////////////////////////////////////////////
// check scalar version
//////////////////////////////////////////////////

    err1 = 0.0;
    Lam  = lam;

    for (n=0; n<N; n++) {
      poisscinv_bisection_scalar(n, N, Lam, Ulo_h, Uhi_h);

      err1 += 0.5*fabs( (Ulo_ex[n]-Ulo_h[n]) + (Uhi_ex[n]-Uhi_h[n]));

      if (n>1) {
        if ( Uhi_h[n-1] <= Ulo_ex[n] ) {
          printf("\n error: lam = %f, n = %d, Uhi_h[n-1] = %20.16g, Ulo_ex[n] = %20.16g \n",
                            lam,      n,      Uhi_h[n-1],           Ulo_ex[n]);
          exit(1);
        }
        if ( Uhi_ex[n-1] <= Ulo_h[n] ) {
          printf("\n error: lam = %f, n = %d, Uhi_ex[n-1] = %20.16g, Ulo_h[n] = %20.16g \n",
		            lam,      n,      Uhi_ex[n-1],           Ulo_h[n]);
          exit(1);
        }
      }
    }

    printf("%10.4g     %9.3g \n",lam,err1);
  }

// free memory 

  free(Ulo_h);
  free(Uhi_h);

  free(ulo_ex);
  free(uhi_ex);
  free(Ulo_ex);
  free(Uhi_ex);
}
