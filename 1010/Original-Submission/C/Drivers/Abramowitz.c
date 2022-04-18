#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<complex.h>
#define Sqr(x) ((x)*(x))
void solve_cubic_analytic(double *coeff, complex double sol[3])
{
  /* solve the cubic coeff[3]*x^3 + coeff[2]*x^2 +  coeff[1]*x + coeff[0] = 0
   * according to the method described in Numerical Recipe book */  
  double a, b, c, Q, R, theta, Q3, R2, A, B;
  const double sqrt32=sqrt(3)/2.0;
  a = coeff[2]/coeff[3];
  b = coeff[1]/coeff[3];
  c = coeff[0]/coeff[3];
  Q = (Sqr(a) - 3.0*b)/9.0;
  R = (2.0*Sqr(a)*a - 9.0*a*b + 27.0*c)/54.0;

  Q3 = Sqr(Q)*Q;
  R2 = Sqr(R);
  if (R2 < Q3)
    {
      theta = acos(R/sqrt(Q3));
      sol[0] = -2.0*sqrt(Q)*cos(theta/3.0)- a/3.0;
      sol[1] = -2.0*sqrt(Q)*cos((theta+2.0*M_PI)/3.0) - a/3.0;
      sol[2] = -2.0*sqrt(Q)*cos((theta-2.0*M_PI)/3.0) - a/3.0;
    }
  else
    {
      A = -copysign(1.0,R)*pow(fabs(R) + sqrt(R2 - Q3),1.0/3.0);
      if (A==0.0)
	B=0.0;
      else
	B = Q/A;
      sol[0] = (A+B) - a/3.0;
      sol[1] = -0.5*(A+B)-a/3.0+I*sqrt32*(A-B);
      sol[2] = -0.5*(A+B)-a/3.0-I*sqrt32*(A-B);
    }
}

void csolve_quartic_abramovitz_cmplx(double *coeff, complex double sol[4])
{
  /* analytic solution of the quartic equation 
   * coeff[4]*x^4 + coeff[3]*x^3 + coeff[2]*x^2 + coeff[1]*x + coeff[0] = 0
   * according to Abramowitz and Stegun book */
  double a3, a2, a1, a0, a32, a12;
  double cb[4], y1;
  double complex solc[3], R, D, E, A, B;
  a3 = coeff[3]/coeff[4];
  a2 = coeff[2]/coeff[4];
  a1 = coeff[1]/coeff[4];
  a0 = coeff[0]/coeff[4];
  a32 = Sqr(a3);
  a12 = Sqr(a1);
  cb[3] = 1.0;
  cb[2] = -a2;
  cb[1] = a1*a3-4.0*a0;
  cb[0] = 4.0*a2*a0-a12-a32*a0;
  solve_cubic_analytic(cb, solc);
  y1 = creal(solc[0]);
  R = csqrt(0.25*a32-a2+y1);
  if (R==0)
    {
      A = 0.75*a32-2.0*a2;
      B = 2.0*csqrt(Sqr(y1)-4*a0);
    }
  else
    {
      A = 0.75*a32-R*R-2.0*a2;
      B = 0.25*(4.0*a3*a2-8.0*a1-a32*a3)/R;
    }  
  D = csqrt(A+B);
  E = csqrt(A-B); 
  sol[0] = -0.25*a3 + 0.5*R + 0.5*D;
  sol[1] = -0.25*a3 + 0.5*R - 0.5*D;
  sol[2] = -0.25*a3 - 0.5*R + 0.5*E;
  sol[3] = -0.25*a3 - 0.5*R - 0.5*E;
}
