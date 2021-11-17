#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include <complex.h>
#define Sqr(x) ((x)*(x))
#ifndef CMPLX
#define CMPLX(x,y) ((x)+I*(y))
#endif
void solve_cubic_analytic_depressed_handle_inf(double b, double c, double *sol)
{
  double Q, R, theta, A, B;
  double QR, QRSQ, KK, sqrtQ, RQ;
  const double PI2=M_PI/2.0, TWOPI=2.0*M_PI;
  Q = -b/3.0;
  R = 0.5*c;
  if (R==0)
    {
      if (b <= 0)
	{
	  *sol=sqrt(-b);
	}
      else
	{
	  *sol=0;
	}
      return;
    }
  
  if (fabs(Q) < fabs(R))
    {
      QR=Q/R;
      QRSQ=QR*QR; 
      KK=1.0 - Q*QRSQ;
    }
  else
    {
      RQ = R/Q;
      KK = copysign(1.0,Q)*(RQ*RQ/Q-1.0);
    }

  if (KK < 0.0)
    {
      sqrtQ=sqrt(Q);
      theta = acos((R/fabs(Q))/sqrtQ);
      if (theta < PI2) 
	*sol = -2.0*sqrtQ*cos(theta/3.0);
      else 
	*sol = -2.0*sqrtQ*cos((theta+TWOPI)/3.0);
    }
  else
    {
      if (fabs(Q) < fabs(R))
	A = -copysign(1.0,R)*cbrt(fabs(R)*(1.0+sqrt(KK)));
      else
	{
	  A = -copysign(1.0,R)*cbrt(fabs(R)+sqrt(fabs(Q))*fabs(Q)*sqrt(KK));
	}
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      *sol = A+B;
    }
}
void solve_cubic_analytic_depressed(double b, double c, double *sol)
{
  double Q, R, theta, Q3, R2, A, B, sqrtQ;
  Q = -b/3.0;
  R = 0.5*c;
  if (fabs(Q) > 1E102 || fabs(R) > 1E154)
    {
      solve_cubic_analytic_depressed_handle_inf(b, c, sol);
      return;
    }
  Q3 = Sqr(Q)*Q;
  R2 = Sqr(R);
  if (R2 < Q3)
    {
      theta = acos(R/sqrt(Q3));
      sqrtQ=-2.0*sqrt(Q);
      if (theta < M_PI/2) 
	*sol = sqrtQ*cos(theta/3.0);
      else 
	*sol = sqrtQ*cos((theta+2.0*M_PI)/3.0);
    }
  else
    {
      A = -copysign(1.0,R)*pow(fabs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      *sol = A+B; /* this is always largest root even if A=B */
    }
}

void calculate_ys(double a, double b, double c, double d, double *phi0)
{
  double g,h,gg,hh;
  double diskr, aq,bq,cq,dq,s;
  int iter;
  double MACHEPS=2.2204460492503131E-16;
  double x, xsq, gx, xxx, maxtt, f, fold, df, xold;
  diskr=9*a*a-24*b;                    
  if(diskr > 0.0)
    { 
      diskr=sqrt(diskr);
      if(a > 0.0)
	s=-2*b/(3*a+diskr);           
      else
	s=-2*b/(3*a-diskr);          
    }
  else
    {      
      s=-a/4;                       
    }
  // !--------------------------- the shift transformation (Horner forms):
  aq=a+4*s;               
  bq=b+3*s*(a+2*s);                             
  cq=c+s*(2*b+s*(3*a+4*s));                      
  dq=d+s*(c+s*(b+s*(a+s)));                      
  gg=bq*bq/9;
  hh=aq*cq;      
  g=hh-4*dq-3*gg;                                
  h=(8*dq+hh-2*gg)*bq/3-cq*cq-dq*aq*aq;         
  solve_cubic_analytic_depressed(g, h, phi0);
  x = *phi0;
  xsq=x*x;
  xxx=x*xsq;
  gx=g*x;
  f = x*(xsq + g) + h;
  if (fabs(xxx) > fabs(gx))
    maxtt = fabs(xxx);
  else
    maxtt = fabs(gx);
  if (fabs(h) > maxtt)
    maxtt = fabs(h);

  if (fabs(f) > MACHEPS*maxtt)
    {
      for (iter=0; iter < 8; iter++)
	{   
	  df =  3.0*xsq + g;
	  if (df==0)
	    {
	      break;
	    }
	  xold = x;
	  x += -f/df;
	  fold = f;
	  xsq = x*x;
	  f = x*(xsq + g) + h;
	  if (f==0)
	    {
	      break;
	    } 
	  if (fabs(f) >= fabs(fold))
	    {
	      x = xold;
	      break;
	    }
    	}
    }
  *phi0 = x;
}

void solve_quadratic_real_coeff(double a, double b, complex double roots[2])
{ 
  double sqrtd, diskr,div,zmax,zmin;
  diskr=a*a-4*b;   
  if(diskr>=0.0)
    {
      if(a>=0.0)
	div=-a-sqrt(diskr);
      else
	div=-a+sqrt(diskr);

      zmax=div/2;

      if(zmax==0.0)
	zmin=0.0;
      else
	zmin=b/zmax;

      roots[0]=CMPLX(zmax,0.0);
      roots[1]=CMPLX(zmin,0.0);
    } 
  else
    {   
      sqrtd = sqrt(-diskr);
      roots[0]=CMPLX(-a/2,sqrtd/2);
      roots[1]=CMPLX(-a/2,-sqrtd/2);      
    }   
}
void solve_quadratic_cmplx_coeff(complex double acx, complex double bcx, complex double sol[2])
{ 
  complex double zx1, zx2, cdiskr;
  cdiskr=csqrt(acx*acx-4.0*bcx);
  zx1 = -0.5*(acx+cdiskr);
  zx2 = -0.5*(acx-cdiskr);
  if (cabs(zx1) > cabs(zx2))
    sol[0] = zx1;
  else
    sol[0] = zx2;
  sol[1]= bcx/sol[0];
} 
void csolve_quartic_shmakov(double *coeff, complex double sol[4])
{
  complex double g1, g2, h1, h2, sq[2];
  double ac, bc, ys, a, b, c, d;

  a=coeff[3]/coeff[4];
  b=coeff[2]/coeff[4];
  c=coeff[1]/coeff[4];
  d=coeff[0]/coeff[4];
  calculate_ys(a, b, c, d, &ys);
  ac = -a;
  bc = (2.0/3.0)*b - ys;
  solve_quadratic_real_coeff(ac, bc, sq);
  g1=sq[0];
  g2=sq[1];
  ac = -(ys+b/3.0);
  bc = d;
  solve_quadratic_real_coeff(ac, bc, sq);
  h1=sq[0];
  h2=sq[1];
  if (cimag(g1)==0 && cimag(h1)==0)
    {
      solve_quadratic_real_coeff(creal(g1), creal(h1), sq);
      sol[0] = sq[0];
      sol[1] = sq[1];
    }
  else
    {
      solve_quadratic_cmplx_coeff(g1, h1, sq);
      sol[0] = sq[0];
      sol[1] = sq[1];
    }
  if (cimag(g2)==0 && cimag(h2)==0)
    {
      solve_quadratic_real_coeff(creal(g2), creal(h2), sq);
      sol[2] = sq[0];
      sol[3] = sq[1];
    }
  else
    {
      solve_quadratic_cmplx_coeff(g2, h2, sq);
      sol[2] = sq[0];
      sol[3] = sq[1];
    }
}

