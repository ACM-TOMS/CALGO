/**********************************************************************************

  DELTAGAMMAINC Fast and Accurate Evaluation of a Generalized Incomplete Gamma
  Function. Copyright (C) 2016 Remy Abergel (remy.abergel AT gmail.com), Lionel
  Moisan (Lionel.Moisan AT parisdescartes.fr).

  This file is a part of the DELTAGAMMAINC software, dedicated to the
  computation of a generalized incomplete gammafunction. See the Companion paper
  for a complete description of the algorithm.

  ``Fast and accurate evaluation of a generalized incomplete gamma function''
  (Rémy Abergel, Lionel Moisan), preprint MAP5 nº2016-14, revision 1.

  This program is free software: you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************/

/* Definition of the G-function
 * ----------------------------
 *
 * We define the function G : (p,x) --> R as follows
 *
 * if x <= p:
 *
 * G(p,x) = exp(x-p*log(|x|)) * integral of s^{p-1} * exp(-sign(x)*s) ds from s = 0 to |x|
 *
 * otherwise:
 *
 * G(p,x) = exp(x-p*log(|x|)) * integral of s^{p-1} * exp(-s) ds from s = x to infinity
 *
 * where
 *
 *   + p is a positive real number
 *   + x is a real number, eventually equal to +infinity.
 *
 *
 * Definition of the generalized incomplete gamma function
 * -------------------------------------------------------
 *
 * We call generalized incomplete gamma function, and we note I_{x,y}^{mu,p},
 * the integral defined by
 *
 *   I_{x,y}^{mu,p} = integral of s^{p-1} * exp(-mu*s) ds
 *
 * where
 *
 *  + mu is a real number non equal to zero (in general we take mu = 1 or -1 but
 *  any nonzero real number is allowed)
 *
 *  + x and y are two nonnegative real numbers such as 0 <= x <= y <= +infinity,
 *    the setting y=+infinity is allowed only when mu > 0
 *
 *  + p is positive real number, p must be an integer when mu < 0.
 *
 *
 * * Brief description of several modules
 * --------------------------------------
 *
 *   + Gammaln approximates the log of the complete gamma function
 *
 *      Gamma(p) = integral of s^{p-1} e^{-s} ds from s=0 to +infinity.
 *
 *     using Pugh's method (a refinement of Lanczos algorithm).
 *
 *   + G_cfrac_lower evaluates the G-function in the domain x <= p using a
 *     continued fraction.
 *
 *   + G_ibp evaluates the G-function in the domain (x < 0 and |x| < max(1,p-1))
 *     using a recursive integration by parts relation (WARNING: this function
 *     cannot be used when mu > 0).
 *
 *   + G_cfrac_upper evaluates the G-function in the domain x > p using a
 *     continued fraction.
 *
 *
 * References
 * ----------
 *
 * + R. Abergel and L. Moisan. 2016. Fast and accurate evaluation of a
 *   generalized incomplete gamma function, preprint MAP5 nº2016-14, revision 1
 *
 * + F. W. J. Olver, D. W. Lozier, R. F. Boisvert, and C. W. Clark
 *   (Eds.). 2010. NIST Handbook of Mathematical Functions. Cambridge University
 *   Press. (see online version at [[http://dlmf.nist.gov/]])
 *
 * + W. H. Press, S. A. Teukolsky, W. T. Vetterling, and
 *   B. P. Flannery. 1992. Numerical recipes in C: the art of scientific
 *   computing (2nd ed.).
 *
 * + G. R. Pugh, 2004. An analysis of then Lanczos Gamma approximation (phd
 *   thesis)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernel.h"

/* plim: compute plim(x), the frontier of the partition of the domain (p,x)
 * detailed in our paper.
 *
 *            |      x              if   0 < x
 *            |
 * plim(x) = <       0              if -9 <= x <= 0
 *            |
 *            | 5.*sqrt(|x|)-5.     otherwise
 * */
static double plim(double x)
{
  return (x >= 0.) ? x : ((x >= -9.) ? 0. : 5.*sqrt(-x)-5.);
}

/* *
 * Gammaln: compute log(Gamma(p)) for p > 0,
 *
 * where Gamma(p) = integral of s^{p-1} e^{-s} ds from s=0 to +infinity.
 *
 * The evaluation is done using the Pugh's method (approximation with 11 terms),
 * which is a refinement of the Lanczos method (approximation with 6 terms).
 *
 * input : p is a positive real number.
 * output : set *gln = log(Gamma(p))
 * */
void Gammaln(double *gln, double p)
{
  double sum,z=p-1.;
  double d[11]={
    2.48574089138753565546e-5, 1.05142378581721974210e+0,
    -3.45687097222016235469e+0, 4.51227709466894823700e+0,
    -2.98285225323576655721e+0, 1.05639711577126713077e+0,
    -1.95428773191645869583e-1, 1.70970543404441224307e-2,
    -5.71926117404305781283e-4, 4.63399473359905636708e-6,
    -2.71994908488607703910e-9 };
  int k;

  for(sum=d[0],k=1; k<=10; k++) sum += d[k]/(z+(double)k);
  *gln = log(1.860382734205265717*sum) -(z+0.5) + (z+0.5L)*log(z+11.400511);
}


/* *
 * G_cfrac_lower: compute G(p,x) in the domain x <= p using a continued fraction
 *
 * inputs :
 *
 *   + p: positive real number
 *   + x: a real number such as x <= p
 *
 * output : set *Gcfrac = G(p,x)
 * */
void G_cfrac_lower(double *Gcfrac, double p, double x)
{
  double c,d,del,f,an,bn;
  int k,n;

  /* deal with special cases */
  if(x==0) { *Gcfrac = 0; return; }

  /* evaluate the continued fraction using Modified Lentz's method. However, as
   * detailed in our paper, we perform manually the first pass (n=1), of the
   * initial Modified Lentz's method. */
  an=1.; bn=p; // set an = a{1}, bn = b{1} (b{1} is non-zero)
  f=an/bn; c=an/DPMIN; d=1./bn; n=2;
  do {
    k=n/2; an = ((n&1) ? k : -(p-1+k))*x; bn++;
    d=an*d+bn; if (d == 0) d=DPMIN;
    c=bn+an/c; if (c == 0) c=DPMIN;
    d=1.0/d; del=d*c; f*=del; n++;
  }
  while((fabs(del-1.0) >= EPS) && (n < ITMAX));
  *Gcfrac = f;

  return;
}


/* *
 * G_ibp: compute G(p,x) in the domain (x<0) using a recursive integration by
 * parts
 *
 * inputs :
 *
 *   + p : positive *integer*
 *   + x : a negative real number such as |x| < max(1,p-1)
 *
 * output : set *Gibp = G(p,x)
 * */
void G_ibp(double *Gibp, double p, double x)
{
  double gln,t,tt,c,d,s,del;
  int l;
  char odd,stop;

  /* initialization */
  Gammaln(&gln,p); t = fabs(x); tt = 1./(t*t);
  odd = ((int)p)%2;

  /* main loop */
  c = 1./t; d = (p-1.); s = c*(t-d); l = 0;
  do {
    c *= d*(d-1)*tt;
    d -= 2.;
    del = c*(t-d);
    s += del; l++;
    stop = (fabs(del) < fabs(s)*EPS);
  }
  while((l<floor((p-2)/2)) && !stop);
  if(odd && !stop) { s += d*c/t; }
  *Gibp = (((odd)?-1.:1.)*exp(-t+gln-(p-1)*log(t))+s)/t;

  return;
}


/* *
 * G_cfrac_upper: compute G(p,x) in the domain x > p using a continued fraction
 *
 * inputs :
 *
 *   + p: positive real number
 *   + x: a positive number x > p, eventually x=infinity
 *
 * output : set *Gcfrac = G(p,x)
 * */
void G_cfrac_upper(double *Gcfrac, double p, double x)
{
  double c,d,del,f,an,bn;
  int i,n;
  char t;

  /* deal with special cases */
  if(isinf(x)) { *Gcfrac = 0; return; }

  /* evaluate the continued fraction using Modified Lentz's method. However, as
   * detailed in our paper, we perform manually the first pass (n=1), of the
   * initial Modified Lentz's method. */
  an=1.,bn=x+1.0-p; // set an = a{1}, bn = b{1}
  t = (bn != 0);
  if(t) { // b{1} is non-zero
    f=an/bn; c=an/DPMIN; d=1./bn; n=2;
  }
  else { // b{1}=0 but b{2} is non-zero, we compute Mcfrac = a{1}/f with f = a{2}/(b{2}+) a{3}/(b{3}+) ...
    an=-(1-p); bn=x+3.0-p;
    f=an/bn; c=an/DPMIN; d=1./bn; n=3;
  }
  i = n-1;
  do {
    an = -i*(i-p); bn += 2.0;
    d=an*d+bn; if (d == 0) d=DPMIN;
    c=bn+an/c; if (c == 0) c=DPMIN;
    d=1.0/d; del=d*c; f*=del; i++; n++;
  }
  while((fabs(del-1.0) >= EPS) && (n < ITMAX));
  *Gcfrac = (t) ? f : 1./f;

  return;
}


/* *
 * G_func: compute G(p,x) using the appropriate routine according to the value
 * of (p,x).
 *
 * inputs :
 *
 *  + p: a positive real number
 *  + x: a real number (eventually x=+inf)
 *
 * output : set *G = G(p,x)
 * */
void G_func(double *G, double p, double x)
{
  if (p >= plim(x)) G_cfrac_lower(G,p,x);
  else if (x < 0) G_ibp(G,p,x);
  else G_cfrac_upper(G,p,x);
}


/* perform 1 iteration of the Romberg approximation of I_{x,y}^{mu,p} */
static void romberg_iterations(R,sigma,n,x,y,mu,p,h,pow2)
     double *R,sigma,x,y,mu,p,h,pow2;
     int n;
{
  int j,m,adr0,adr0_prev;
  double sum,xx,pow4;

  adr0_prev = ((n-1)*n)/2;
  adr0 = (n*(n+1))/2;

  /* initialization */
  for(sum=0,j=1;j<=pow2;j++) {
    xx = x + ((y-x)*(2*j-1))/(2.*pow2);
    sum += exp(-mu*xx+(p-1)*log(xx)-sigma);
  }
  R[adr0] = 0.5*R[adr0_prev] + h*sum;

  /* main loop */
  pow4=4.;
  for(m=1;m<=n;m++) {
    R[adr0+m] = (pow4*R[adr0+(m-1)]-R[adr0_prev+(m-1)])/(pow4-1.);
    pow4 *= 4.;
  }

  return;
}

/* Estimate I_{x,y}^{mu,p} using a Romberg approximation: this algorithm
 * computes rho and sigma such as the Romberg approximation of I_{x,y}^{mu,p}
 * is given by I_{x,y}^{mu,p} = to rho * exp(sigma)                    */
void romberg_estimate(rho,sigma,x,y,mu,p)
     double *rho,*sigma,x,y,mu,p;
{
  int n,adr0;
  double *R,h,pow2,relerr,relneeded;
  char cont;

  /* memory allocation */
  R = (double*) malloc (((NITERMAX_ROMBERG+1)*(NITERMAX_ROMBERG+2))/2*sizeof(double));

  /* initialization (n=1) */
  *sigma = -mu*y + (p-1)*log(y);
  R[0] = 0.5*(y-x)*(exp(-mu*x+(p-1)*log(x)-(*sigma))+1.);

  /* loop for n > 0 */
  relneeded = EPS/TOL_ROMBERG;
  adr0 = 0;
  n = 1;
  h = (y-x)/2.; // n=1, h = (y-x)/2^n
  pow2 = 1;     // n=1; pow2 = 2^(n-1)
  if (NITERMAX_ROMBERG >=1)
    do {
      // update R //
      romberg_iterations(R,*sigma,n,x,y,mu,p,h,pow2);
      h /= 2.;
      pow2 *= 2;
      // estimate relative error //
      adr0 = (n*(n+1))/2;
      relerr = fabs((R[adr0+n]-R[adr0+n-1])/R[adr0+n]);
      // check for stopping criterion //
      n++; cont = ((n<=NITERMAX_ROMBERG) && (relerr > relneeded));
    } while (cont);

  /* save output (Romberg estimate) and free memory */
  *rho = R[adr0+(n-1)];
  free(R);

  return;
}


/* deltagammainc: our algorithm for the approximation of I_{x,y}^{mu,p}.
 *
 * inputs :
 *
 *  + mu is a real number different from 0
 *
 *  + x and y are two nonnegative real numbers such as 0 <= x <= y <= +inf,
 *    and the setting y=+infinity is allowed only when mu > 0
 *
 *  + p is positive real number, p must be an integer when mu < 0
 *
 * outputs: this procedure computes (rho,sigma,method) described below:
 *
 * (rho,sigma) are such as the approximated value of I_{x,y}^{mu,p} is
 * I=(*rho)*exp(*sigma)
 *
 * (*method) is a flag describing how I was computed: (*method)==1 when I is
 * estimated using a difference, (*method)==2 when I is estimated using
 * Romberg's method, (*method) < 0 in other (trivial) cases (like I=0, or I =
 * Gamma(p)).
 * */
void deltagammainc(double *rho, double *sigma, char *method, double x, double y, double mu, double p)
{
  double mA,mB,mx,my,nA,nB,nx,ny,gln;

  /* deal with particular cases */
  if(x==y){ *rho = 0.; *sigma = -INFINITY; (*method) = -1; return; }
  if (x==0. && isinf(y)){ *rho = 1.; Gammaln(sigma,p); (*sigma) = (*sigma)-p*log(mu); (*method) = -2; return; }

  /* initialization */
  G_func(&mx,p,mu*x); nx = (isinf(x)) ? -INFINITY : -mu*x + p*log(x);
  G_func(&my,p,mu*y); ny = (isinf(y)) ? -INFINITY : -mu*y + p*log(y);

  /* KERNEL: Compute (mA,nA) and (mB,nB) such as I_{x,y}^{mu,p} can be
     approximated by the difference A-B, where A >= B >= 0, A = mA*exp(nA) and B
     = mB*exp(nB). When the difference involves more than one digit loss due to
     cancellation errors, the integral I_{x,y}^{mu,p} is evaluated using the
     Romberg approximation method. */
  if(mu<0) { mA = my; nA = ny; mB = mx; nB = nx; }
  else {
    if(p<plim(mu*x)) { mA = mx; nA = nx; mB = my; nB = ny; }
    else if (p < plim(mu*y)) {
      Gammaln(&gln,p);
      mA = 1.; nA = gln-p*log(mu);
      nB = fmax(nx,ny);
      mB = mx*exp(nx-nB) + my*exp(ny-nB);
    }
    else { mA = my; nA = ny; mB = mx; nB = nx; }
  }
  /* compute (rho,sigma) such as rho*exp(sigma) = A-B */
  *rho = mA-mB*exp(nB-nA);
  *sigma = nA;

  /* if the difference involved a significant loss of precision, compute Romberg estimate */
  if(!isinf(y) && ((*rho)/mA < TOL_DIFF)) {
    romberg_estimate(rho,sigma,x,y,mu,p);
    (*method) = 2;
  }
  else (*method) = 1;
}


/* convert (rho,sigma) representation into scientific notation, that is, set
 * (a,b) such as rho*exp(sigma) = a * 10^(b) where a in [0,10) and b has an
 * integer value
 * */
void scientific_notation(double *a, double *b, double rho, double sigma)
{
  double l10,delta;

  l10 = log(10);
  delta = log10(rho)+sigma/l10;
  *b = floor(delta);
  *a = exp((delta-(*b))*l10);

}
