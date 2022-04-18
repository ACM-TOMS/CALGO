#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kernel.h"

extern void G_func(double*,double,double); // see source file 'kernel.c'
extern void deltagammainc(double*,double*,char*,double,double,double,double); // see source file 'kernel.c'

/* tabulated values of (rho,sigma) such as I(x,y,mu,p) = rho*exp(sigma), for several values of (x,y,mu,p) */
static double p[15]  = { 1, 5,10,12,14,  1,  5, 10, 20, 1, 3,10, 1,10,20};
static double x[15]  = { 9, 9, 9, 9, 9,100,100,100,100, 5, 5, 5,20,20,20};
static double y[15]  = {11,11,11,11,11,120,120,120,120,10,10,10,25,25,25};
static double mu[15] = { 1, 1, 1, 1, 1,  1,  1,  1,  1,-1,-1,-1,-1,-1,-1};
static double Iref[15] = {1.06708103296433890e-04, 9.56616980230235669e-01, 8.95942017652358167e+04, 8.93104948155385003e+06, 9.02034141170810289e+08, 3.72007596835318789e-44, 3.87343326443145790e-36, 4.08366058817001986e-26, 4.57980828029277503e-06, 2.18780526357041399e+04, 1.80364717146940695e+06, 1.12951155494984621e+13, 7.15197341419760822e+10, 2.00168223708455573e+23, 1.47339480836645219e+37};
static double rel_Iref[15] = {6.17e-16, 1.39e-15, 1.30e-14, 2.79e-15, 2.45e-15, 4.88e-15, 3.90e-15, 7.89e-15, 1.58e-14, 1.08e-15, 1.45e-16, 3.88e-15, 1.61e-15, 8.64e-16, 1.49e-15};

/* tabulated values of G(p,x) for several values of (p,x) (see below) */
static double Gref[25] = {4.33274822761946129e-01, 1.09773581433105088e-02, 1.00907250703918066e-03, 1.00090072050430255e-04, 1.00009000720050403e-05, 1.10975931173828892e-02, 1.28772193213532231e-01, 1.10974291120600913e-03, 1.00999797021514453e-04, 1.00099097114199831e-05, 1.01009071471357288e-03, 1.11097413973301546e-03, 3.99699388464565804e-02, 1.11097397244931695e-04, 1.01009070421084097e-05, 1.00100090071046924e-04, 1.01009997959810302e-04, 1.11109739572180804e-04, 1.25665794460617313e-02, 1.11109739404560889e-05, 1.00010000900071005e-05, 1.00100099097104160e-05, 1.01010090704106744e-05, 1.11110973938932034e-05, 3.96666393667642864e-03};
static double rel_Gref[25] = {5.31e-16, 9.11e-18, 9.10e-17, 5.12e-17, 1.60e-16, 1.56e-16, 2.12e-15, 1.02e-16, 8.26e-17, 5.93e-17, 1.07e-16, 6.91e-16, 3.16e-15, 1.90e-16, 4.22e-17, 4.48e-17, 8.60e-17, 1.96e-16, 1.46e-14, 1.63e-16, 2.53e-16, 2.34e-18, 2.28e-17, 4.28e-18, 5.47e-14};

int main(int argc, char **argv)
{
  int id;
  double rho,sigma,rel,rel_ref,pp,xx,G,I;
  char method,status,STATUS=1;

  /* Unit tests for I */
  printf("\nUnit tests for module 'deltagammainc': compute (rho,sigma) such as\n");
  printf("rho * exp(sigma) = I = integral over [x,y] of s^{p-1} * exp(-mu*s) ds.\n\n");
  printf("(Check for results displayed in Table V of the companion paper).\n\n");
  for(id=0;id<15;id++) {
      deltagammainc(&rho,&sigma,&method,x[id],y[id],mu[id],p[id]);
      I = rho*exp(sigma);
      //rel = fabsl(1.L-(long double)I/(long double)Iref[id]); // use long double to reduce effect of cancellation in the computation of the relative error
      rel = fabs(1-I/Iref[id]);
      rel_ref = rel_Iref[id];
      status = (rel<=10*rel_ref);
      STATUS &= status;
      printf("+ (p=%-2.2g, x=%-3.3g, y=%-3.3g, mu=%-2.2g) : I=%-23.17e, relative error (expected~%0.0e, actual=%0.0e), status = %s\n",p[id],x[id],y[id],mu[id],I,rel_ref,rel,status?"PASSED":"FAILED");
  }
  printf("\n");

  printf("\nUnit tests for module 'Gfunc': compute G(p,x) such as\n");
  printf("if x <= p: G(p,x) = exp(x-p*log(|x|)) * integral over [0,|x|] of s^{p-1} * exp(-sign(x)*s) ds\n");
  printf("otherwise: G(p,x) =  exp(x-p*log(x))  * integral over [x,inf] of s^{p-1} * exp(-s) ds.\n\n");

  /* Unit tests for G */
  for(id=0,pp=10;pp<=1e5;pp*=10){
    for(xx=10;xx<=1e5;xx*=10,id++) {
      G_func(&G,pp,xx);
      //rel = fabs(1.L-(long double)G/(long double)Gref[id]); // same comment as above
      rel = fabs(1-G/Gref[id]);
      rel_ref = rel_Gref[id];
      rel = fmax(1.22e-16,rel);
      rel_ref = fmax(1.22e-16,rel_ref);
      status = (rel<=10*rel_ref);
      STATUS &= status;
      printf("+ (p=%0.0e, x=%0.0e) : G=%23.17e, relative error (expected~%0.0e, actual=%0.0e), status = %s\n",pp,xx,G,rel_ref,rel,status?"PASSED":"FAILED");
    }
  }
  printf("\n");

  /* Final message */
  if(STATUS) printf("All unit tests were successfully passed!\n\n");
  else printf("Some of the unit tests failed!\n\n");

  return EXIT_SUCCESS;
}
