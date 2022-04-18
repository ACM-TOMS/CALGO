/* Algorithm 409: Discrete Chebychev Curve Fit, H. Schmitt,
  Communications of the ACM, May 1971. Translated to C, April 2020.
  Copyright © 1971 Association for Computing Machinery, Inc. */
#include <math.h>
#include <stdbool.h>
#define entier (int)floor
static int sign(double x) {return (x>0)-(x<0);}

int acm409(int m, int n, int k, // Renamed from ‘approx’.
    double const x[], double const y[],
    double epsh, int * maxit, int ref[],
    double * hmax, double h[], double a[])
/* This procedure computes the best approximation poly-
  nomial in the sense of Chebychev of required degree m to a set
  of n distinct points given by their abscissas and ordinates (array
  x, y [1:n]).  The abscissas must be arranged in increasing order
  x[1] < x[2] < ··· < x[n]. The desired polynomial is even, odd,
  or mixed for k = 2, k = 1, or k = 0, respectively. It is expected
  that x[1] ≥ 0 in case of k = 2 and x[1] > 0 in case of k = 1.
  Leveling according to the exchange method described by Stiefel
  [1] is done up to a tolerance of abs(epsh). The sign of epsh
  decides whether ref is expected to supply entry data (cf. param-
  eter ref).
    maxit enters an upper limit for the number of exchange steps
  allowed and returns the number of steps actually performed.
  The parameter ref is used to carry entry data only if epsh < 0. It
  is an integer array containing the subscripts of the points to be
  used as initial reference. The lower array bound is 1, the upper
  bound (say p) is m + 2 in the case of mixed (k = 0) polynomials,
  entier ((m+3)/2) in the case of odd (k = 1), and entier
  ((m+4)/2) in the case of even (k = 2) polynomials. It is expected
  that 1 ≤ ref[1] < ref[2] < ··· < ref[p] ≤ n. Unless an initial
  reference is not explicitly given by means of the array ref and
  indicated by a value epsh < 0, the points lying next to the so-
  called Chebychev abscissas (with regard to the interval [x[1],
  x[n]]) are determined to start off the algorithm. As output, this
  parameter returns the reference belonging to the approximation
  polynomial.
    The output parameters are hmax to return the maximum devia-
  tion, an array h[1:n] to return the approximation errors at all
  given points, and an array a[0:m] to carry the polynomial
  coefficients. The array h containing the approximation errors is
  introduced as a formal parameter to allow a drawing of the error
  function to be made outside the procedure. This provides a means
  to look at the quality of the computed approximation and is
  recommended to the user. A totally leveled approximation
  polynomial should have an error function with well charac-
  terized extrema of equal height.
    Three emergency exits are provided for extraordinary events.
  exparameter is an exit to be used when entry data are entered
  incorrectly, exmaxit is used when the best fit is not found within
  the maximum number of exchange steps allowed. In this case,
  the parameter ref denotes a new reference which may be used as
  entry data for a further call of approx. The exit exsign is used
  when the approximation errors at the points of reference do not
  alternate in sign. In this case, accuracy of the computer is insuffi-
  cient to generate an approximation polynomial of the required
  degree.
    Acknowledgment. The author wishes to express his apprecia-
  tion to Prof. Dr. W. Barth for many valuable discussions on the
  subject of Chebychev approximation.
Reference
  1. Stiefel, E. L. Numerical methods of Chebychev approximation.
  In On Numerical Approximation, R. Langer, (Ed.), U. Wis-
  consin Press, 1958, pp. 217–232 */
{
  int i, j, p, q1, q2, r;  bool k0, k1;
  k0 = k == 0;  k1 = k == 1;
  q1 = k1? 1 : 0;
  q2 = k0? 1 : 2;
  for (i = 0; i <= m; ++i) a[i] = 0;
  if (! k0) m = entier((m-q1) * 0.5 + 0.1);
  p = m + 2;
  // Check for properly given parameters:
  if (n < p || m < 0 || (! k0 && (! k1 || x[1] <= 0)
    && (k != 2 || x[1] < 0))) goto exparameter;
  for (i = 2; i <= n; ++i)
    if (x[i] <= x[i-1]) goto exparameter;
  {
#ifdef __cplusplus
    auto exchange=[](int n, int p, double h[], double epsh, int z[])
#else
    bool exchange(int n, int p, double h[], double epsh, int z[])
#endif
    /* This procedure performs the exchange technique.
      The number of points and the number of reference points
      are entered by n and p. The approximation errors at different
      points are compared relative to epsh. The subscripts of the
      points of reference are carried by z[1] ··· z[p] of the integer
      array z[0:p+1], a parameter which serves to enter the
      former and return the new reference, z[0] and z[p+1] are
      for internal use only and are expected to have the values 0
      and n + 1. If both the old and new references are equal to
      each other, a jump to the label equal occurs. No global
      quantities are contained within this procedure. */
    {
      int i, j, l, index=0, indl=0, indr=0, sig, ze;
      double hz1, hzp, max, maxl, maxr;
      l = 0;  sig = -sign(h[z[1]]);
      if (sig == 0) sig = 1;
      for (i = 1; i <= p; ++i)
      {
        max = 0;  sig = -sig;  ze = z[i+1] - 1;
        for (j = z[i-1] + 1; j <= ze; ++j)
        if ((h[j]-max) * sig > 0)
        { max = h[j];  index = j; }
        if (fabs (max-h[z[i]]) > fabs(max) * epsh)
        { z[i] = index;  l = 1; }
      }
      maxl = maxr = 0;
      for (j = z[p] + 1; j <= n; ++j)
      if ((maxr -h[j]) * sig > 0)
      { maxr = h[j];  indr = j; }
      hz1 = h[z[1]];  sig = sign(hz1);
      for (j = 1; j <= z[1] - 1; ++j)
      if ((maxl-h[j]) * sig > 0)
      { maxl = h[j];  indl = j; }
      maxl = fabs(maxl);  maxr = fabs(maxr);
      hz1 = fabs(hz1);  hzp = fabs(h[z[p]]);
      if (l == 0)
      {
        if (maxl - hzp <= maxl * epsh &&
        maxr - hz1 <= maxr * epsh) goto equal;
      }
      if (maxl == 0 && maxr == 0) goto end;
      if (maxl > maxr)
      {
        if (maxl > hzp) goto shl;
        else if (maxr >= hz1) goto shr;
      }
      else
      {
        if (maxr > hz1) goto shr;
        else if (maxl >= hzp) goto shl;
      }
      goto end;
shr:
      index = z[1];
      for (i = 1; i <= p - 1; ++i) z[i] = z[i+1];
      z[p] = indr;
      if (maxl > 0)
      for (i = 1; i <= p - 1; ++i)
      {
        if (fabs (h[indl]) >= fabs (h[z[i]]))
        { j = z[i];  z[i] = indl;  indl = index;
          index = j; }
        else goto end;
      }
      goto end;
shl:
      index = z[p];
      for (i = p; i >= 2; --i) z[i] = z[i-1];
      z[1] = indl;
      if (maxr > 0)
      for (i = p; i >= 2; --i)
      {
        if (fabs (h[indr]) >= fabs(h[z[i]]))
        { j = z[i];  z[i] = indr;  indr = index;
          index = j; }
        else goto end;
      }
end:  return false;  equal:return true;
    }; // procedure exchange
    double arg, max, pi, q, s, t, dt, x1 , xa, xe;  bool b1, b2;
    double xx[1+n], aa[1+m], daa[1+m], c[1+p], d [1+p];
    int z[1+p+1];
    // Set up of initial reference:
    z[0] = 0;  z[p+1] = n + 1;
    if (epsh < 0)
    {
      j = 0;
      for (i = 1; i <= p; ++i)
      {
        r = z[i] = ref[i];
        if (j < r) j = r; else goto exparameter;
      }
      if (j > n) goto exparameter;
      epsh = fabs (epsh);  goto m1;
    }
    pi = 3.14159265;  x1 = x[1];  xe = x[n];
    if (k0)
    { xa = xe + x1;  xe = xe - x1;
      arg = pi/(m+1); }
    else
    { xa = 0;  xe = xe + xe;
      arg = pi/(2*(m+1)+q1); }
    for (j = p; j >= 1; --j)
    {
      x1 = xa + xe * cos (arg * (p-j));  r = z[j+1];
      for (i = r-1; i >= 2; --i)
      if (x[i] + x[i-1] <= x1) goto m0;
      i = 1;
m0:
      z[j] = r > i? i : r - 1;
    }
    if (z[1] >= 1) goto m1;
    for (j = 1; z[j] < j; j = j + 1) z[j] = j;
m1:
    for (i = 0; i <= m; ++i) aa[i] = 0;
    for (i = 1; i <= n; ++i)
    { h[i] = y[i];  q = x[i];
      xx[i] = k0? q : q * q;
    }
    b1 = b2 = false;  r = -1;  t = 0;
iterat:
    r = r + 1;  s = 1.0;
    // Computation of the divided difference schemes:
    if (k1)
    {
      for (i = 1; i <= p; ++i)
      {
        s = -s;  j = z[i];  q = x[j];
        c[i] = (h[j] + s * t)/q;  d[i] = s/q;
      }
    }
    else
    for (i = 1; i <= p; ++i)
    { s = -s;  c[i] = h[z[i]] + s * t;  d[i] = s; }
    for (i = 2; i <= p; ++i)
      for (j = p; j >= i; --j)
    {
      q = xx[z[j]] - xx[z[1+j-i]];
      c[j] = (c[j] - c[j-1])/q;
      d[j] = (d[j] - d[j-1])/q;
    }
    dt = -c[p]/d[p];  t = t + dt;
    // Computation of the polynomial coefficients:
    for (i = m; i >= 0; --i)
    {
      daa[i] = c[i+1] + dt * d[i+1];  q = xx[z[i+1]];
      for (j = i; j <= m - 1; ++j)
        daa[j] = daa[j] - q * daa[j+1];
    }
    for (i = 0; i <= m; ++i) aa[i] = aa[i] + daa[i];
    // Evaluation of the polynomial to get the approxima-
    //  tion errors:
    max = 0;
    for (i = 1; i <= n; ++i)
    {
      s = aa[m];  q = xx[i];
      for (j = m - 1; j >= 0; --j) s = s * q + aa[j];
      if (k1) s = s * x[i];
      q = h[i] = y[i] - s;
      if (fabs (q) > max) max = fabs(q);
    }
    // Test for alternating signs:
    j = -sign (h[z[1]]);
    for (i = 2; i <= p; ++i)
      if (sign (h[z[i]]) == j) j = -j; else
      { b1 = true;  goto m2; }
    // Search for another reference:
    if (exchange (n, p, h, epsh, z)) goto m2;
    if (r < *maxit) goto iterat; else b2 = true;
    // Results to output parameters:
m2:
    for (i = 0; i <= m; ++i) a[q1+i*q2] = aa[i];
    for (i = 1; i <= p; ++i) ref[i] = z[i];
    *hmax = max;  *maxit = r;
    if (b1) goto exsign;
    if (b2) goto exmaxit;
  } return 0;  exmaxit:return 1;  exsign:return 2;  exparameter:return 3;
} // procedure acm409
