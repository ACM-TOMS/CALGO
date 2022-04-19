/*
   --------------------------------------------------------------
   Testexample  of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
*/


#include "adouble.h"
#include "adutils.h"
#include <stream.h>
#include <stdio.h>

#ifdef __GNUG__
#include <std.h>
#else
#include <stdlib.h>
#endif

/* 
   This program computes first order directional derivatives 
   for the helmholtz energy function.
*/
    struct tms *buffer;
    int n,deg,ideg;
    adouble y,z,u;
    adouble den,r,t0,t1,t2,t3;
    adouble ttim[6];
    adouble ttim2[6];
    adouble bx,te,he,xax,tem;
 
adouble energy(int n,adouble x[], adouble bv[])
{
  adouble r,he;
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



void main( int argc,char *argv[]) 
{
  int nf,n,j;
  if (argc < 2)
    {
      printf("Error- I expect at least one integer: #of independents/10 \n");
      exit(-1);
    }
  nf = atoi(argv[1]); 
  double r;
  double result=0.0;
  n = 10*nf; 
  adouble* x=new adouble[n+1];
  adouble* bv=new adouble[n+1];
  adouble* dir=new adouble[n+1];
  double* grad=new double[n+1];
  r=  1.0/n;
  for (j=0; j < n; j++) 
    {
      dir[j]=0.1*j;
      bv[j]= 0.02*(1.0+fabs(sin(double(j))));
    }
  
  int imd_rev =1;
  trace_on(1,imd_rev);
  for (j=0; j < n; j++) 
    {
      x[j]<<=r*sqrt(1.0+j);
    }
  he = energy(n,x,bv);
  he >>= result;

  trace_off();
  printf("%f -- energy\n",result);
  reverse(1,1,n,0,1.0,grad);
  for (int l=0; l < n;l++)
    {
      printf("%d, %f,  \n",l,grad[l]);
    }
  printf("%f -- energy\n",result);
  free((char*) x);
  free((char*) bv);
  free((char*) dir);
  free((char*) grad);
}

