#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<complex.h>
extern void oqs_quartic_solver_cmplx(complex double coeff[5], complex double sol[4]);
extern void oqs_quartic_solver(double coeff[5], complex double sol[4]);
int main(void)
{
  double c[5]; 
  complex double cc[5], roots[4], xc[4];
  int k;
  xc[0] = -1+3.0*I;
  xc[1] = 1-I;
  xc[2] = -2.5+3*I;
  xc[3] = 1-2*I;
  cc[4] = 1.0;
  cc[3] = -(xc[0]+xc[1]+xc[2]+xc[3]);
  cc[2] = xc[0]*xc[1] + (xc[0]+xc[1])*(xc[2]+xc[3]) + xc[2]*xc[3]; 
  cc[1] = -xc[0]*xc[1]*(xc[2]+xc[3]) - xc[2]*xc[3]*(xc[0]+xc[1]);
  cc[0] = xc[0]*xc[1]*xc[2]*xc[3];

  printf("Complex Quartic:\n");
  printf("x^4+(%.16G+(%.16G)*I)*x^3+(%.16G+(%.16G)*I)*x^2+(%.16G+(%.16G)*I)*x+(%.16G+(%.16G)*I)==0\n",
	 creal(cc[3]),cimag(cc[3]),creal(cc[2]),cimag(cc[2]),creal(cc[1]),cimag(cc[1]),
	 creal(cc[0]),cimag(cc[0]));
  printf("\nExact Roots:\n");
  for (k=0; k < 4; k++)
    {
      printf("root #%d: %.16G+(%.16G)*I\n", k, creal(xc[k]),cimag(xc[k]));
    }
  printf("\nCalculated Roots:\n");
  oqs_quartic_solver_cmplx(cc,roots);
  for (k=0; k < 4; k++)
    {
      printf("root #%d: %.16G+(%.16G)*I\n", k, creal(roots[k]),cimag(roots[k]));
    }

  xc[0] = 1E4+3*I;
  xc[1] = 1E4-3*I;
  xc[2] = -7+1E3*I;
  xc[3] = -7-1E3*I;
  c[4] = 1.0;
  c[3] = -creal(xc[0]+xc[1]+xc[2]+xc[3]);
  c[2] = creal(xc[0]*xc[1] + (xc[0]+xc[1])*(xc[2]+xc[3]) + xc[2]*xc[3]); 
  c[1] = creal(-xc[0]*xc[1]*(xc[2]+xc[3]) - xc[2]*xc[3]*(xc[0]+xc[1]));
  c[0] = creal(xc[0]*xc[1]*xc[2]*xc[3]);

  printf("\n=========================================================================== \n");
  printf("\nReal Quartic:\n");
  printf("x^4+(%.16G)*x^3+(%.16G)*x^2+(%.16G)*x+(%.16G)==0\n", c[3],c[2],c[1],c[0]);
  printf("\nExact Roots:\n");
  for (k=0; k < 4; k++)
    {
      printf("root #%d: %.16G+(%.16G)*I\n", k, creal(xc[k]),cimag(xc[k]));
    }
  printf("\nCalculated Roots:\n");
  oqs_quartic_solver(c,roots);
  for (k=0; k < 4; k++)
    {
      printf("root #%d: %.16G+(%.16G)*I\n", k, creal(roots[k]),cimag(roots[k]));
    }
  return 0;
}
