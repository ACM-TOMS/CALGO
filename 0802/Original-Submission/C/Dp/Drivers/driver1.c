/* ALC2DLIB-Version 1.0 implemented by WH 22.3.99, all rights reserved */
/*please report problems or bugs to whoer@statistik.wu-wien.ac.at      */
/*                               or hormannw@boun.edu.tr               */

/* This file:               main1.c

   Demonstrates, how to call functions of  "alc2d.c" to generate
   random pairs of two-dimensional log-concave distributions

   To compile this example on a unix-computer just type

   cc alc2d.c main1.c -lm

*/

#include <math.h>
#include <stdio.h>
#include "alc2d.h"


/****************************************************************************/
/* Include your random number generator; the generator here is only
   suitable for small simulations!!!!  */

static unsigned long int urn=0;/*seed of the uniform generator*/

double unif()/*uniform r n Marsaglia m=2^32 a=69069 c=1*/
/*uniform random number generator*/
{
   urn=(69069*(urn)+1);
   return(urn/4.294967296e9);
}
/*******************************************************************/
/*code the logarithm of the density function you want to sapmle from:
  as this is the first example we take the standard normal distribution;
  allthough we do not use parameters, we must icnlude them in the prototype
  as the generation algorithm is expecting them as well*/

/*log of bivariate standard normal density*/
double hstandardnormal(double x[2],double par[])
{return -0.5*(x[0]*x[0]+x[1]*x[1]);}

/*partial derivative with respect to x[0]*/
double hnx(double x[2],double par[])
{return(-x[0]);}

/*partial derivative with respect to x[1]*/
double hny(double x[2],double par[])
{return(-x[1]);}

/*******************************************************************/
int main()
{ long int i; 
  double x[2],sp[1][2],eq[4][3];
  void *t1=NULL;

/*first we have to choose a starting-value close to the mode, e.g.:*/
  sp[0][0]=0.1; sp[0][1]=0.2;

/* As the normal distribution is defined on the whole R^2, we need not
   define the domain of the distribution,
   but as we have chosen only one starting-value we have to use
   a bounded auxiliary domain: In this example we use the square
   with the corners (-1,-1), (-1,1), (1,-1), and (1,1), as auxiliary
   domain.
   It is bordered by the four lines   0=1+x; 0=1+y; 0=-1+x; 0=-1+y
   lines are stored in the form:
                    0 = eq[0][0] + eq[0][1]*x + eq[0][2]*y ;*/

  eq[0][0]=1.; eq[0][1]=1.;eq[0][2]=0.;
  eq[1][0]=1.; eq[1][1]=0.;eq[1][2]=1.;
  eq[2][0]=-1.; eq[2][1]=1.;eq[2][2]=0.;
  eq[3][0]=-1.; eq[3][1]=0.;eq[3][2]=1.;


  t1=setup(1,20,sp,0,eq,4,hstandardnormal,hnx,hny,NULL,0);
/* 1 ... we have one starting point
   20... the algorithm can use up to 20 design points
   sp... array containing the starting points
   0 ... no of equalities for the domain (as it is the whole R^2)
   eq... array containing the equalities for domain and auxiliary domain
   4 ... no of equalities for the auxiliary domain
   hstandardnormal ... logarithm of the density function
   hnx... partial derivative with respect to x
   hny... partial derivative with respect to y
   NULL.. no array of parameters is needed
   0  ... no of parameters of the distribution*/


/* now we can sample from the distribution, the result is stored in x[]*/
  for(i=0;i<1000;i++) sample2d(x,t1);

  printf("The last generated pair was the pair (%f|%f)\n",x[0],x[1]);

/* after generation the allocated memory must be freed*/
  freesetup(t1);

  exit(0);
}

