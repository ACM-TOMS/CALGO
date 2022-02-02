/*
   --------------------------------------------------------------
   Testexample  of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
*/

#include <sys/types.h>
#include <stdio.h>
#include <sys/times.h>
#include <math.h>
#include <stdlib.h>
#define delta 0.000001

/* 
   This program computes first order directional derivatives 
   for the helmholtz energy function.
*/

    struct tms *buffer;
    int n,deg,ideg;
    double y,z,u;
    double den,r,t0,t1,t2,t3;
    double ttim[6];
    double ttim2[6];
    double bx,te,he,xax,tem;
 
double energy(int n,double x[],double bv[])
{
  double r,he;
  int i,j;
  xax = 0;
  bx = 0;
  he =0;
  te =0;
  
  for (i=0; i < n; i++)
    {
      he += x[i]*log(x[i]);
      bx +=  bv[i]*x[i];
      tem = (2.0/(1.0+i+i))*x[i];
      for (j=0; j<i; j++) 
	tem += (1.0/(1.0+i+j))*x[j];
      xax += x[i]*tem;
    };
  xax *= 0.5;
  he = 1.3625E-3*(he-te*log(1.0-bx));
  r=sqrt(2.0);
  he = he - log((1+bx*(1+r))/(1+bx*(1-r)))*xax/bx;
  
  return he;
}


void main(int argc,char *argv[]) 
{
  int nf;
  int n,j,l;
  double r;
  double he2;
  double he1;
  double q;
  double j_prime;
  double *x;
  double *bv;
  double *dir;

  if (argc < 2)
    {
      printf("Error- I expect at least one integer: #of independents/10 \n");
      exit(-1);
    }
  nf = atoi(argv[1]);
  n = 10*nf;
  x  = (double *)malloc((n+1)*sizeof(double));
  bv = (double *)malloc((n+1)*sizeof(double));
  dir = (double *)malloc((n+1)*sizeof(double));
  r=1.0/n;
  for (j=0; j < n; j++) 
    {
      dir[j]=0.1*j;
      j_prime = j;
      bv[j]= 0.02*(1.0+fabs( sin( j_prime ) ) );
      /* for (int i=0; i < n; i++) 
	   am[i][j] = 1.0/(1.0+i+j); */
    }
  for (j=0; j < n; j++) 
    {
      x[j]=r*sqrt(1.0+j);
    }
  he = energy(n,x,bv); 
  he2 = he;
  printf("%f -- energy\n",he);
  for (l=0; l < n;l++)
    {
      
      x[l]=x[l]+delta;
      he1 = energy(n,x,bv);
      x[l]=x[l]-delta;
      q = ((he1-he)/delta);
      printf("%d, %f,  \n",l,q);
    }
  free(dir);
  free(bv);
  free(x);
}

