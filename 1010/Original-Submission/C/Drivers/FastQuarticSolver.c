#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include <complex.h>
#include <float.h>
#define USE_ABRAMOWITZ
#define NBUF 12
#define Sqr(x) ((x)*(x))
extern void csolve_quartic_abramovitz_cmplx(double *coeff, complex double sol[4]);
extern void solve_quadratic(double coeff[3], int *numsol, double *sol);
extern void solve_quadratic_cmplx(double coeff[3], complex double *sol);
extern void solve_cubic_analytic(double *coeff, complex double sol[3]);
void solve_quadratic(double coeff[3], int *numsol, double *sol)
{
  /* numeric error safe version of solve_quadratic from Numerical Recipe */
  double delta, a, b, c, q;
  a = coeff[2];
  b = coeff[1];
  c = coeff[0];
  delta = Sqr(b) - 4.0*a*c;
  if (delta > 0.0)
    {
      q = -0.5*(b+copysign(1.0,b)*sqrt(delta));
      sol[0] = q/a;
      sol[1] = c/q;
      *numsol = 2;
    } 
  else if (delta == 0)
    {
      sol[0] = -b/(2.0*a);
      *numsol = 1;
    }
  else
    {
      *numsol = 0;
    }
}

int fqs_cmp_func(const void* aa, const void *bb)
{
  complex double ai, bi;
  ai = *((complex double *)aa);
  bi = *((complex double *)bb);
  if (cabs(ai) > cabs(bi))
    return -1;
  else if (cabs(ai)==cabs(bi))
    return 0;
  else
    return 1;
}
void csolve_quartic_ferrari_cmplx(double *coeff, complex double sol[4])
{
  /* questa soluzione di fatto è quella di abramovitz */
  double a43, a2a3a4, a44, a4, a3, a2, a1, a0, a32, a42, a3a4;
  double cb[4], p, q, r, m=0, cq[3];
  double complex solc[3], solq[2];
  complex double sm, A, B, C, Dp, Dm;
  const double sqrt2=1.4142135623730950488016887242097; 
  int k, mzero;
  a4 = coeff[4];
  a3 = coeff[3];
  a2 = coeff[2];
  a1 = coeff[1];
  a0 = coeff[0];
  if (a3==0 && a2==0 && a1==0 && a0 ==0)
    {
      sol[0] = 0.0;
      sol[1] = 0.0;
      sol[2] = 0.0;
      sol[3] = 0.0;
      return;
    }
  a32 = Sqr(a3);
  a42 = Sqr(a4);
  a3a4=-0.25*a3/a4;
  a43 = a42*a4;
  a2a3a4 = a2*a3*a4;
  a44 = a42*a42;
  p = (8.0*a2*a4-3.0*a32)/8.0/a42;
  q = (a32*a3 - 4.0*a2a3a4 + 8.0*a1*a42)/8.0/(a43);
  r = (-3.0*a32*a32+256.0*a0*a43-64.0*a1*a3*a42+16*a3*a2a3a4)/256.0/(a44);
  cb[3] = 8.0;
  cb[2] = 8.0*p;
  cb[1] = 2*Sqr(p)-8.0*r;
  cb[0] = -Sqr(q);
  solve_cubic_analytic(cb, solc);
  mzero=1;
  for (k=0; k < 3; k++)
    {
      if (cimag(solc[k])==0)
	{
	  if (creal(solc[k])!=0)
	    m = creal(solc[k]);
	  mzero=0;
	  break;
	}
    }
  if (mzero)
    {
      /* hence q=0 and quartic is a biquadratic */
      cq[2] = 1.0;
      cq[1] = p;
      cq[0] = r;
      solve_quadratic_cmplx(cq, solq);
      sol[0] = csqrt(solq[0]);
      sol[1] = -csqrt(solq[0]);
      sol[2] = csqrt(solq[1]);
      sol[3] = -csqrt(solq[1]);
      return;
    }
  sm = csqrt(m);
  A = 0.5*sqrt2*sm;
  B = sqrt2*q/sm;
  C = 2.0*(p+m);
  Dp = 0.5*csqrt(-(C + B));
  Dm = 0.5*csqrt(-(C - B));
  sol[0] = a3a4 + A + Dp;
  sol[1] = a3a4 + A - Dp;
  sol[2] = a3a4 - A + Dm;
  sol[3] = a3a4 - A - Dm;
}
int eps_identical(double *eps)
{
  int j;
  for (j=1; j < NBUF; j++)
    {
      /* if actual eps is identical to one of the four previous ones then terminate! */
      if (eps[0]==eps[j])
	return 1;
    }
  // if it gets here they are all equal!
  return 0;
}
double nummax=0, totcall=0;
int fqsits=0;
void backward_optimizer(double *alpha, double *beta, double *gamma, double *delta, double a, double b, double c, double d, int *kchosen)
{
  double e1[2], e2[2], e3[2], e4[2];
  double U23[2], U33[2], L43[2], U44[2], x1[2], x2[2], x3[2], x4[2], y1[2], y2[2], y3[2], y4[2];
  double eps[2][NBUF];
  const int MAXITS=17; 
  int ignore[2];
  int k, j, its;
  totcall++;
  for (k=0; k < NBUF; k++)
    {
      eps[0][k] = eps[1][k] = 0;
    }
  ignore[0]=ignore[1]=0;
  for (k=0; k < 2; k++)
    {
      e1[k] = a - alpha[k] - gamma[k];
      e2[k] = b - beta[k] - alpha[k]*gamma[k] - delta[k];
      e3[k] = c - beta[k]*gamma[k] - alpha[k]*delta[k];
      e4[k] = d - beta[k]*delta[k];
    }
  for (its=0; its < MAXITS; its++)
    {
      for (k=0; k < 2; k++)
	{
	  if (ignore[k])
	    continue;
	  U23[k] = alpha[k] - gamma[k];
	  U33[k] = beta[k] - delta[k] - gamma[k]*U23[k];
	  L43[k] = -delta[k]*U23[k]/U33[k];
	  U44[k] = beta[k] - delta[k] - L43[k]*U23[k];
	  x1[k] = e1[k];
	  x2[k] = e2[k] - gamma[k]*x1[k];
	  x3[k] = e3[k] - delta[k]*x1[k] - gamma[k]*x2[k];
	  x4[k] = e4[k] - delta[k]*x2[k] - L43[k]*x3[k];
	  if (U44[k]==0.0|| U33[k]==0.0)
	    {
	      ignore[k]=1;
	      if (ignore[0]==1 && ignore[1]==1)
		{
		  *kchosen = k;
		  return;
		}
	    }
	  y4[k] = x4[k]/U44[k];
	  y3[k] = (x3[k]-U23[k]*y4[k])/U33[k];
	  y2[k] = x2[k] - U23[k]*y3[k] - y4[k];
	  y1[k] = x1[k] - y3[k];
	  alpha[k] = alpha[k] + y1[k];
	  beta[k] = beta[k] + y2[k];
	  gamma[k] = gamma[k] + y3[k];
	  delta[k] = delta[k] + y4[k];
	  e1[k] = a - alpha[k] - gamma[k];
	  e2[k] = b - beta[k] - alpha[k]*gamma[k] - delta[k];
	  e3[k] = c - beta[k]*gamma[k] - alpha[k]*delta[k];
	  e4[k] = d - beta[k]*delta[k];
	  /* shift epsilon's */ 
	  for (j=NBUF-1; j > 0; j--)
	    eps[k][j] = eps[k][j-1];
	  eps[k][0] = fabs(e1[k])+fabs(e2[k])+fabs(e3[k])+fabs(e3[k]);
	  // convergence
	  if (eps[k][0] == 0.0)
	    {
	      *kchosen=k;
	      fqsits = its;
	      return;
	    }
	  // cyclic condition
	  else if (eps_identical(eps[k]))
	    {
	      fqsits = its;
	      *kchosen=k;
	      return;
	    }
	}
    }
  fqsits = its;
  nummax++;
  if (ignore[0]==1)
   {
     *kchosen=1;
     return;
   }
  if (ignore[1]==1)
    {
      *kchosen = 0;
      return;
    }
      
  if (eps[0][0] < eps[1][0])
    {
      *kchosen = 0;
    }
  else
    {
      *kchosen = 1;
    }
}
int error_handler1(double a, double b, double c, double d, complex double solqua[4])
{
  double cq[3], alpha, beta, eps1, eps2;
  complex double sq[2];
  int k;
  alpha = a*0.5;
  beta = (b - Sqr(alpha))*0.5;
  eps1 = c - 2.0*(alpha)*(beta);
  eps2 = d - Sqr(beta);
  if (eps1==0 && eps2==0)
    {
      /* now we solve the quadratic equation providing the roots of the original quartic */
      cq[2] = 1.0;
      cq[1] = alpha;
      cq[0] = beta;
      solve_quadratic_cmplx(cq, sq);
      for (k=0; k < 2; k++)
	{
	  solqua[k] = solqua[k+2] = sq[k];
	}
      return 1;
    }
  return 0; 
}
int error_handler2(double a, double b , double c, double d, int *numsol, complex double solqua[4])
{
  double cq[3], x1[2], x2[2], eps1, eps2;
  int nsq, k;
  cq[2] = 1.0;
  cq[1] = a*0.5;
  cq[0] = b/6.0;
  solve_quadratic(cq, &nsq, x1);
  if (nsq==0)
    {
      *numsol=0;
      return 0;
    } 
  x2[0] = -a - 3.0*x1[0];
  x2[1] = -a - 3.0*x1[1];
     
  for (k=0; k < 2; k++)
    {
      eps1 = c + Sqr(x1[k])*(x1[k]+3.0*x2[k]);
      eps2 = d - Sqr(x1[k])*x1[k]*x2[k];
      if (eps1==0 && eps2==0)
	{
	  *numsol=2;
	  solqua[0] = solqua[1] = solqua[2] = x1[k];
	  solqua[3] = x2[k];	  
	  return 1;
	}
      
    }
  return 0;
   
}

void initial_guess_fast_quart_solver(double *alpha, double *beta, double *gamma, double *delta, double a, double b, double c, double d, complex double csol[4])
{
  double coeff[5];
  double phi1, phi2, c1, c2, L1, L2, L3, y1, y2;
  int k;
  coeff[4] = 1.0;
  coeff[3] = a;
  coeff[2] = b;
  coeff[1] = c;
  coeff[0] = d;
#ifdef USE_ABRAMOWITZ
  csolve_quartic_abramovitz_cmplx(coeff, csol);
#elif defined(USE_FERRARI)
  csolve_quartic_ferrari_cmplx(coeff, csol);
#elif defined(USE_SALZER) 
  solve_fourth_deg_cmplx(coeff, csol);
#endif
  qsort(csol, 4, sizeof(complex double), fqs_cmp_func);
  alpha[0] = -creal(csol[0]+csol[1]);
  beta[0] = creal(csol[0]*csol[1]);
  alpha[1] = -creal(csol[1]+csol[2]);
  beta[1] = creal(csol[1]*csol[2]);
  for (k=0; k < 2; k++)
    {
      phi1 = 1.0 + Sqr(alpha[k])+Sqr(beta[k]);
      phi2 = alpha[k]*(1.0+beta[k]);
      c1 = a - alpha[k] + alpha[k]*(b-beta[k])+beta[k]*c;
      c2 = b - beta[k] + alpha[k]*c + beta[k]*d;
      L1 = sqrt(phi1);
      L3 = phi2/L1;
      L2 = sqrt(phi1-phi2*phi2/phi1);
      y1 = c1/L1;
      y2 = (c2 - y1*L3)/L2;
      delta[k] = y2/L2;
      gamma[k] = (y1 - delta[k]*L3)/L1;
    }
}
void fast_quartic_solver(double coeff[5], complex double solqua[4])
{
  /* solution of the quartic equation 
   * coeff[4]*x^4 + coeff[3]*x^3 + coeff[2]*x^2 + coeff[1]*x + coeff[0] = 0
   * according to the fast quartic solver proposed in Ref. [38] */
  
  double alpha[2], beta[2], gamma[2], delta[2];
  double a, b, c, d;
  double cq[3];
  complex double csol[4], csq[2];
  int k, setchosen, numsol; 
  a = coeff[3]/coeff[4];
  b = coeff[2]/coeff[4];
  c = coeff[1]/coeff[4];
  d = coeff[0]/coeff[4];

  initial_guess_fast_quart_solver(alpha, beta, gamma, delta, a, b, c, d, csol);

  fqsits = 0;
  if (error_handler1(a, b, c, d, solqua))
    return;

  if (error_handler2(a, b, c, d, &numsol, solqua))
    return;
  backward_optimizer(alpha, beta, gamma, delta, a, b, c, d, &setchosen);
  /* now we solve the two quadratic equation providing the four roots of the original quartic */
  cq[2] = 1.0;
  cq[1] = alpha[setchosen];
  cq[0] = beta[setchosen];
  solve_quadratic_cmplx(cq, csq);
  for (k=0; k < 2; k++)
    {
      solqua[k] = csq[k];
    }
  cq[2] = 1.0;
  cq[1] = gamma[setchosen];
  cq[0] = delta[setchosen];
  solve_quadratic_cmplx(cq, csq);
  for (k=0; k < 2; k++)
    {
      solqua[k+2] = csq[k];
    }
}
