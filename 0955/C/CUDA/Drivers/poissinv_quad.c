#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <quadmath.h>


//////////////////////////////////////////////////
// compute reference solution in quad precision
//////////////////////////////////////////////////

void poissinv_quad(int N, float lam, float  *u_lo, float  *u_hi,
                                     double *U_lo, double *U_hi) {
  int n;
  __float128 q_lam, q_s, q_n; 

  char buf[128];

  q_lam = lam;
  q_s   = 0.0q;

  for (n=0; n<N; n++) {
    // doing it this way avoids problems with limited range
    q_n  = (__float128) n;
    q_s += expq( -q_lam + q_n*logq(q_lam) - lgammaq(q_n+1.0q) );
    if (q_s > 1.0q - 1.0e-1000q) q_s = 1.0q - 1.0e-1000q;

//
// single precision interval bisection
//
    u_hi[n] = 1.0f;
    u_lo[n] = 0.0f;
    float u_mid = 0.5f*(u_hi[n] + u_lo[n]);
 
    while (u_mid>u_lo[n] & u_mid<u_hi[n]) {
      if ((__float128) u_mid > q_s)
        u_hi[n] = u_mid;
      else
        u_lo[n] = u_mid;

      u_mid = 0.5f*(u_hi[n] + u_lo[n]);
    }
 
//
// double precision interval bisection
//
    U_hi[n] = 1.0;
    U_lo[n] = 0.0;
    double U_mid = 0.5*(U_hi[n] + U_lo[n]);
 
    while (U_mid>U_lo[n] & U_mid<U_hi[n]) {
      if ((__float128) U_mid > q_s)
        U_hi[n] = U_mid;
      else
        U_lo[n] = U_mid;

      U_mid = 0.5*(U_hi[n] + U_lo[n]);
    }
 
  }
}


void poisscinv_quad(int N, float lam, double *U_lo, double *U_hi) {
  int n;
  __float128 q_lam, q_s, q_n; 

  char buf[128];

  q_lam = lam;
  q_s   = 0.0q;

  for (n=N-1; n>=0; n--) {
    // doing it this way avoids problems with limited range
    q_n  = (__float128) n;
    q_s += expq( -q_lam + q_n*logq(q_lam) - lgammaq(q_n+1.0q) );
    if (q_s > 1.0q - 1.0e-1000q) q_s = 1.0q - 1.0e-1000q;

//
// double precision interval bisection
//
    U_hi[n] = 1.0;
    U_lo[n] = 0.0;
    double U_mid = 0.5*(U_hi[n] + U_lo[n]);
 
    while (U_mid>U_lo[n] & U_mid<U_hi[n]) {
      if ((__float128) U_mid > q_s)
        U_hi[n] = U_mid;
      else
        U_lo[n] = U_mid;

      U_mid = 0.5*(U_hi[n] + U_lo[n]);
    }
 
  }
}
