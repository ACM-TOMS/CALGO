/* ALC2DLIB-Version 1.0 implemented by WH 22.3.99, all rights reserved */
/*please report problems or bugs to whoer@statistik.wu-wien.ac.at      */
/*                               or hormannw@boun.edu.tr               */

/* This file:               main2.c

   Demonstrates, how to call functions of  "alc2d.c" to generate
   random pairs of several two-dimensional log-concave distributions

   To compile this example on a unix-computer just type

   cc alc2d.c main2.c -lm

*/

#include <math.h>
#include <stdio.h>
#include "alc2d.h"



#define MIN(a,b) ((a)<(b)?a:b)

/****************************************************************************/
static unsigned long int urn=0;/*seed of the uniform generator*/

double unif()/*uniform r n Marsaglia m=2^32 a=69069 c=1*/
/*uniform random number generator*/
{
   urn=(69069*(urn)+1);
   return(urn/4.294967296e9);
}

/*******************************************************************/
/*Two-dimensional normal distribution*/
/*par[0]..R, par[1]..S1, par[2]..S2*/

double hn(double x[2],double par[]) /*log of bivariate normal density*/
{return
(-0.5/(1.-par[0]*par[0])*(x[0]*x[0]*(1./(par[1]*par[1]))-
    2.*par[0]*x[0]*x[1]*(1./(par[1]*par[2]))+x[1]*x[1]*(1./(par[2]*par[2]))));}

double hnx(double x[2],double par[]) /*partial derivative w.r.t. x*/
{return(-0.5/(1.-par[0]*par[0])*
	(2*x[0]*(1./(par[1]*par[1]))-2.*par[0]*x[1]*(1./(par[1]*par[2]))));}

double hny(double x[2],double par[]) /*partial derivative w.r.t. x*/
{ return(-0.5/(1.-par[0]*par[0])*
	 (-2.*par[0]*x[0]*(1./(par[1]*par[2]))+2*x[1]*(1./(par[2]*par[2]))));}

/*******************************************************************/
/*"cut normal distribution" of the paper*/
/*log of normal distribution "cut through" by a plane
  by taking the minimum of the paraboloid and the plane*/
/*par[0]..R, par[1]..S1, par[2]..S2, par[3]...NCA, par[4]...NCB, par[5]...NCC*/

double hnc(double x[2],double par[]) /*log of bivariate cut normal density*/

{ double fn;
  fn=(-0.5/(1.-par[0]*par[0])*(x[0]*x[0]*(1./(par[1]*par[1]))-
		     2.*par[0]*x[0]*x[1]*(1./(par[1]*par[2]))+
			       x[1]*x[1]*(1./(par[2]*par[2]))));
  return MIN(par[5]+par[3]*x[0]+par[4]*x[1],fn);
}

double hncx(double x[2],double par[]) /*partial derivative w.r.t. x*/
{ double fn,fe;
  fn=(-0.5/(1.-par[0]*par[0])*(x[0]*x[0]*(1./(par[1]*par[1]))-
		     2.*par[0]*x[0]*x[1]*(1./(par[1]*par[2]))+x[1]*x[1]*(1./(par[2]*par[2]))));
  fe=par[5]+par[3]*x[0]+par[4]*x[1];
  if(fn<fe)
    return(-0.5/(1.-par[0]*par[0])*(2*x[0]*(1./(par[1]*par[1]))-
				    2.*par[0]*x[1]*(1./(par[1]*par[2]))));
  else return par[3];
}

double hncy(double x[2],double par[]) /*partial derivative w.r.t. y*/
{ double fn,fe;
  fn=(-0.5/(1.-par[0]*par[0])*(x[0]*x[0]*(1./(par[1]*par[1]))-
		     2.*par[0]*x[0]*x[1]*(1./(par[1]*par[2]))+x[1]*x[1]*(1./(par[2]*par[2]))));
  fe=par[5]+par[3]*x[0]+par[4]*x[1];
  if(fn<fe)
    return(-0.5/(1.-par[0]*par[0])*(-2.*par[0]*x[0]*(1./(par[1]*par[2]))+
				    2*x[1]*(1./(par[2]*par[2]))));
  else return par[4];
}
/*******************************************************************/
/*Distribution NS2 of the paper*/
/*par[0]...n*/

double href(double x[],double par[]) /* log of density*/
{ double r;
  r=sqrt(x[0]*x[0]+x[1]*x[1]);

  return MIN(par[0]*par[0]-par[0]*r,0);
}
double hrefx(double x[],double par[]) /*partial derivative w.r.t. x*/
{ double r;
  r=sqrt(x[0]*x[0]+x[1]*x[1]);

  if(r>(par[0])) return(-par[0]*x[0]/r);
  else return 0.;
}
double hrefy(double x[],double par[]) /*partial derivative w.r.t. y*/
{ double r;
  r=sqrt(x[0]*x[0]+x[1]*x[1]);

  if(r>(par[0])) return(-par[0]*x[1]/r);
  else return 0.;
}
/*******************************************************************/
/*distribution NS3 of the paper
  par[0]..a, par[1]..b, par[2]..c, par[3]..d, par[4]..e,  
  if a,b,d,e>0 and 4bd>=c^2 the distribution is log concave on R^2*/


double hq(double a[2],double par[]) /*log of bivariate quartic density*/
{ double x2,y2,x,y;
  x=a[0];
  y=a[1];
  x2=x*x;
  y2=y*y;
  return(-(par[0]*x2*x2+par[1]*x2+par[2]*x*y+par[3]*y2+par[4]*y2*y2));
}
double hqx(double a[2],double par[]) /*partial derivative w.r.t. x*/
{ double x,y;
  x=a[0];
  y=a[1];
return(-(4*par[0]*x*x*x+2*par[1]*x+par[2]*y));
}

double hqy(double a[2],double par[]) /*partial derivative w.r.t. y*/
{ double x,y;
  x=a[0];
  y=a[1];
  return(-(par[2]*x+2*par[3]*y+4*par[4]*y*y*y));
}
/*******************************************************************/
/*bivariate beta-distribution*/

/*par[0]..a1, par[1]..a2, par[2]..a3*/

double hb(double x[2],double par[])
/*log of density of bivariate beta distribution
  domain of the distribution is the triangle (0,0), (0,1) (1,0)*/
{return
(par[0]-1.)*log(x[0])+(par[1]-1.)*log(x[1])+(par[2]-1.)*log(1.-x[0]-x[1]);}

double hbx(double x[2],double par[]) /*partial derivative w.r.t. x*/
{return (par[0]-1.)/x[0]-(par[2]-1.)/(1.-x[0]-x[1]);}

double hby(double x[2],double par[]) /*Partial derivative w.r.t. y*/
{return (par[1]-1.)/x[1]-(par[2]-1.)/(1.-x[0]-x[1]);}

/*function can be used to initialize the domain of the beta-distribution*/
void betainit(double eq[][3])
{ /*line x=0*/
  eq[0][0]=0.;
  eq[0][1]=1.;
  eq[0][2]=0.;
  /* line y=0 */
  eq[1][0]=0.;
  eq[1][1]=0.;
  eq[1][2]=1.;
  /* line x+y=1 */
  eq[2][0]=-1.;
  eq[2][1]=1.;
  eq[2][2]=1.;
}
/*******************************************************************/
/*bivariate gamma-distribution
  Becker and Roux 1981: South African Statistical Journal 15, 1-12*/
/*Caution!!! This distribution is not log-concave for all choices*/
/*  of the parameters, allthough it is allways logconcave
  for a1, a2>1 for x< y and y<x. But along the line x=y there can
  occur local minima!!*/

/*par[0]..a1, par[1]..a2, par[2]..l1, par[3]..l2, par[4]...b1, par[5]...b2*/

double hgam(double x[2],double par[])
/*log of density of bivariate gamma distribution x[0] and x[1]>=0*/
{
if(x[0]<=x[1])
    return (par[0]-1.)*log(x[0])+(par[1]-1.)*log(par[3]*(x[1]-x[0])+x[0])-
           (1./par[4]+1./par[5]-par[3]/par[5])*x[0]-par[3]*x[1]/par[5];
else
    return (par[1]-1.)*log(x[1])+(par[0]-1.)*log(par[2]*(x[0]-x[1])+x[1])-
           (1./par[4]+1./par[5]-par[2]/par[4])*x[1]-par[2]*x[0]/par[4];

}
double hgamx(double x[2],double par[]) /*partial derivative w.r.t. x*/
{
if(x[0]<=x[1])
    return (par[0]-1.)/x[0]+(par[1]-1.)*(1.-par[3])/(par[3]*(x[1]-x[0])+x[0])
            -(1./par[4]+1./par[5]-par[3]/par[5]);
else 
    return (par[0]-1.)*par[2]/(par[2]*(x[0]-x[1])+x[1])-par[2]/par[4];
}
double hgamy(double x[2],double par[]) /*partial derivative w.r.t. y*/
{
if(x[0]<=x[1])
    return (par[1]-1.)*par[3]/(par[3]*(x[1]-x[0])+x[0])-par[3]/par[5];
else 
    return (par[1]-1.)/x[1]+(par[0]-1.)*(1.-par[2])/(par[2]*(x[0]-x[1])+x[1])
            -(1./par[4]+1./par[5]-par[2]/par[4]);
}


/*function can be used to initialize the domain of the gamma-distribution*/
void gammainit(double eq[][3])
{ /*line x=0*/
  eq[0][0]=0.;
  eq[0][1]=1.;
  eq[0][2]=0.;
  /* line y=0 */
  eq[1][0]=0.;
  eq[1][1]=0.;
  eq[1][2]=1.;
}


/*******************************************************************/
/* distribution NS1 of the paper */

double hg(double x[2],double par[])
/*log of density of a simple 2-dimensional gamma-distribution
  domain is the half-plain with x>0 */
{return log(x[0])-x[0]*x[0]-x[0]*x[1]-x[1]*x[1];}


double hgx(double x[2],double par[]) /*Partielle Ableitung hb nach x*/
{return 1./x[0]-2.*x[0]-x[1];}

double hgy(double x[2],double par[]) /*Partielle Ableitung hb nach y*/
{return -x[0]-2.*x[1];}

/*******************************************************************/
#define HD 4 /*Half length of the side of the square used as auxiliary domain*/

int main()
{ long int i,wid; 
  double x[8][2],sp[1][2],eq[10][3],par[10],pair[2];
  void *t1=NULL,*t2=NULL,*t3=NULL,*t4=NULL,*t5=NULL,*t6=NULL,*t7=NULL,*t8=NULL;


/*auxiliary domain is the square (-HD,HD)x(-HD,HD) */
  eq[0][0]=HD;
  eq[0][1]=1.;
  eq[0][2]=0.;
  eq[1][0]=HD;
  eq[1][1]=0.;
  eq[1][2]=1.;
  eq[2][0]=-HD;
  eq[2][1]=1.;
  eq[2][2]=0.;
  eq[3][0]=-HD;
  eq[3][1]=0.;
  eq[3][2]=1.;


/*We take a point close the mode as starting point*/
  sp[0][0]=0.1;
  sp[0][1]=0.2122;

/*setting the parameters of the normal distribution*/
  par[0]=0.5;/*R*/ 
  par[1]=12.;/*S1*/ 
  par[2]=5.;/*S2*/
  printf("Setup for the normal distribution with r=%f, s1=%f and s2=%f\n",par[0],par[1],par[2]);
  t1=setup(1,50,sp,0,eq,4,hn,hnx,hny,par,3);
/*t1 is now pointer to the structure that holds all information for sampling from the
  normal distribution with R=0.5, S1=12 and S2=5  */

/*As second example we generate from a bivariate normal distribution restricted to the domain
  x>=0 an y>=-x 
  For the definition of the domain we use*/
  eq[0][0]=0.;
  eq[0][1]=1.;
  eq[0][2]=0.;
  eq[1][0]=0.;
  eq[1][1]=1.;
  eq[1][2]=1.;


/*for the auxiliary domain we add the two lines x=3 y=3 */
  eq[2][0]=-3.;
  eq[2][1]=1.;
  eq[2][2]=0.;
  eq[3][0]=-3.;
  eq[3][1]=0.;
  eq[3][2]=1.;


/*We take as starting point*/
  sp[0][0]=1.;
  sp[0][1]=0.;


  par[0]=0.95;/*R*/ 
  par[1]=10.;/*S1*/ 
  par[2]=10.;/*S2*/
  printf("Setup for the normal distribution with r=%f, s1=%f and s2=%f\n",par[0],par[1],par[2]);
  printf("over the domain : x>=0 and y>-x\n");
  t2=setup(1,20,sp,2,eq,4,hn,hnx,hny,par,3);
/* (1,20,sp... shows, that we use one starting-point and that we allow 20 design points;
   ...2,eq,4.. means, that there are 2 equalities for the domain of the distribution
  and 2 additional equalities for an auxiliary domain 
  t2 is now pointer to the structure that holds all information for sampling from the
  normal distribution with R=0.95, S1=10 and S2=10  over the restricted domain*/

/*the NS2 distribution is again generated over the whole R2, so we use the same
  auxiliary domain, as for the unrestricted normal distribution */

/*auxiliary domain is the square (-HD,HD)x(-HD,HD) */
  eq[0][0]=HD;
  eq[0][1]=1.;
  eq[0][2]=0.;
  eq[1][0]=HD;
  eq[1][1]=0.;
  eq[1][2]=1.;
  eq[2][0]=-HD;
  eq[2][1]=1.;
  eq[2][2]=0.;
  eq[3][0]=-HD;
  eq[3][1]=0.;
  eq[3][2]=1.;


/*We take a point close to the mode as starting point*/
  sp[0][0]=0.01;
  sp[0][1]=0.1;


  par[0]=3.;/*n*/
  printf("set-up for distribution NS2 with N=%f\n",par[0]);
  t3=setup(1,50,sp,0,eq,4,href,hrefx,hrefy,par,1);


  par[0]=0.05;/*R*/ par[1]=1.;/*S1*/ par[2]=1.;/*S2*/
  par[3]=0.5;/*A*/ par[4]=1.;/*B*/ par[5]= -1.;/*C*/
  printf("Setup for cut normal-distr.: R=%f S1=%f S2=%f\n",par[0],par[1],par[2]);
  printf("A %f B %f C %f \n",par[3],par[4],par[5]);
  t4=setup(1,50,sp,0,eq,4,hnc,hncx,hncy,par,6);

  par[0]=0.03;/*A*/ par[1]=0.1;/*B*/ par[2]=1.;/*C*/
  par[3]=0.5;/*D*/ par[4]=0.01;/*E*/
/*checks concavity*/
  if(4*par[1]*par[3]<par[2]*par[2])/*case not concave, thus c is changed*/
    par[2]=4*par[1]*par[3];
  printf("\nsetup for NS3-distr.: A=%f B=%f C=%f\n",par[0],par[1],par[2]);
  printf("D %f E %f \n",par[3],par[4]);
  t5=setup(1,50,sp,0,eq,4,hq,hqx,hqy,par,5);


/*to store the domain of the beta distribution*/
  betainit(eq);

  par[0]=20.;/*a1*/ par[1]=3.;/*a2*/ par[2]=3.;/*a3*/

  sp[0][0]=0.87;/* Starting point*/
  sp[0][1]= 0.015;
  printf("setup for beta-distr.: a1=%f a2=%f a3=%f\n",par[0],par[1],par[2]);
  t6=setup(1,50,sp,3,eq,3,hb,hbx,hby,par,3);
/*the ...3,eq,3.. means, that there are 3 equalities for the domain of the distribution
  but no additional equalities for an auxiliary domain */

/*this the domain of the NS-1 distribution*/
 /*line x=0*/
  eq[0][0]=0.;
  eq[0][1]=1.;
  eq[0][2]=0.;


/*auxiliary domain is the rectangle (0,HD)x(-HD,HD) 
  so we have to add the follwing three equalities */
  eq[1][0]=HD;
  eq[1][1]=0.;
  eq[1][2]=1.;
  eq[2][0]=-HD;
  eq[2][1]=1.;
  eq[2][2]=0.;
  eq[3][0]=-HD;
  eq[3][1]=0.;
  eq[3][2]=1.;


  sp[0][0]=0.5;/* Starting point*/
  sp[0][1]= 0.1;
  printf("setup for NS1 (gamma)-distr.: ");
  t7=setup(1,50,sp,1,eq,4,hg,hgx,hgy,NULL,0);
/*the ...1,eq,4.. means, that there is 1 equality for the domain of the distribution
  and 3 additional equalities for an auxiliary domain */
/******************************************************************/
/*to store the domain of the gamma distribution, 2 equalities*/
  gammainit(eq);
/*auxiliary domain is the rectangle (0,20)x(0,20) 
  so we have to add the follwing two equalities */
  eq[2][0]=-10;
  eq[2][1]=1.;
  eq[2][2]=0.;
  eq[3][0]=-10;
  eq[3][1]=0.;
  eq[3][2]=1.;


/*Caution!!! This distribution is not log-concave for all possible parameter
values!!!!!!!!!!!!!!!!!!!!!*/

  par[0]=4./*a1*/; par[1]=3./*a2*/;
  par[2]=3./*l1*/; par[3]=2./*l2*/;
  par[4]=3./*b1*/; par[5]=1./*b2*/;


  sp[0][0]=1.;/* Starting point*/
  sp[0][1]= 1.;
  printf("setup for gamma-distr.: a1=%f a2=%f \n",par[0],par[1]);
  printf("l1=%f l2=%f b1=%f b2=%f\n",par[2],par[3],par[4],par[5]);
  t8=setup(1,50,sp,2,eq,4,hgam,hgamx,hgamy,par,6);
/*the ...2,eq,4.. means, that there are 2 equalities for the domain of the distribution
  and 2 additional equalities for the auxiliary domain */



/* Now we could start a simulation and sample without problems from the eight two-dimensional
   distributions, that were initialized above. */

   printf("We generate 5  matrices , which contain a random pair of each distribution\n");
   printf("The parameters are those that we have specified above\n");

   for(i=0;i<5;i++)
     { sample2d(x[0],t1);
       sample2d(x[1],t2);
       sample2d(x[2],t3);
       sample2d(x[3],t4);
       sample2d(x[4],t5);
       sample2d(x[5],t6);
       sample2d(x[6],t7);
       sample2d(x[7],t8);
       printf("matrix %ld\n",i);
       printf("normal:            (%.2f|%.2f)\n",x[0][0],x[0][1]);
       printf("restricted-dnormal:(%.2f|%.2f)\n",x[1][0],x[1][1]);
       printf("NS2:               (%.2f|%.2f)\n",x[2][0],x[2][1]);
       printf("cut normal:        (%.2f|%.2f)\n",x[3][0],x[3][1]);
       printf("NS3                (%.2f|%.2f)\n",x[4][0],x[4][1]);
       printf("beta:              (%.2f|%.2f)\n",x[5][0],x[5][1]);
       printf("NS1:               (%.2f|%.2f)\n",x[6][0],x[6][1]);
       printf("Gamma:             (%.2f|%.2f)\n",x[7][0],x[7][1]);
     }
  printf("If you do not like the messages of the sample2d-function, that gives information\n");
  printf("about the updated hat, you have to change the define OUTPUT to -1 in alc2d.h\n");

/* No we want to estimate the acceptance probabilities of the different distributions;
   but this is only possible if the define MEMCONTROL is >=1 in alc2d.h */

  wid=50000;

/* First we have to set the repetition counter to 0 */
  setwidcount(0);

  printf("We are generating %ld pairs from the normal distribution\n",wid);
  for(i=0;i<wid;i++) sample2d(pair,t1);
  printf("Observed average acceptance probability d3.1: %f\n",(double)wid/returnwidcount());
  
/*As we will not sample from that distribution any longer, we free the allocated memory*/
  freesetup(t1);

  setwidcount(0);

  printf("We are generating %ld pairs from the restricted normal dist. As we allowed only\n",wid);
  printf("only 20 design points (calling the set-up) we expect a lower acceptance prob.\n");
  for(i=0;i<wid;i++) sample2d(pair,t2);
  printf("Observed average acceptance probability d3.1: %f\n",(double)wid/returnwidcount());
  
  printf("The maximal number of corners of one polygon up to now was %d .\n",returnncmax());

/*As we will not sample from that distribution any longer, we free the allocated memory*/
  freesetup(t2);

  setwidcount(0);

  printf("We are generating %ld pairs from the NS2 distribution\n",wid);
  printf("This function has a big constant part and thus the acceptance-probability is high.\n");
  for(i=0;i<wid;i++) sample2d(pair,t3);
  printf("Observed average acceptance probability d3.1: %f\n",(double)wid/returnwidcount());
  
  printf("The maximal number of corners of one polygon up to now was %d .\n",returnncmax());
  printf("This high number comes from the big, circle-shaped plateau of the NS2-distr.\n");

/*As we will not sample from that distribution any longer, we free the allocated memory*/
  freesetup(t3);

  setwidcount(0);

  printf("We are generating %ld pairs from the cut-normal distribution\n",wid);
  for(i=0;i<wid;i++) sample2d(pair,t4);
  printf("Observed average acceptance probability d3.1: %f\n",(double)wid/returnwidcount());
  
/*As we will not sample from that distribution any longer, we free the allocated memory*/
  freesetup(t4);

  setwidcount(0);

  printf("We are generating %ld pairs from the NS3 distribution\n",wid);
  for(i=0;i<wid;i++) sample2d(pair,t5);
  printf("Observed average acceptance probability d3.1: %f\n",(double)wid/returnwidcount());
  
/*As we will not sample from that distribution any longer, we free the allocated memory*/
  freesetup(t5);

  setwidcount(0);

  printf("We are generating %ld pairs from the beta distribution\n",wid);
  for(i=0;i<wid;i++) sample2d(pair,t6);
  printf("Observed average acceptance probability d3.1: %f\n",(double)wid/returnwidcount());
  
/*As we will not sample from that distribution any longer, we free the allocated memory*/
  freesetup(t6);

  setwidcount(0);

  printf("We are generating %ld pairs from the NS1 distribution\n",wid);
  for(i=0;i<wid;i++) sample2d(pair,t7);
  printf("Observed average acceptance probability d3.1: %f\n",(double)wid/returnwidcount());
  
/*As we will not sample from that distribution any longer, we free the allocated memory*/
  freesetup(t7);

  setwidcount(0);

  printf("We are generating %ld pairs from the Gamma distribution\n",wid);
  for(i=0;i<wid;i++) sample2d(pair,t8);
  printf("Observed average acceptance probability d3.1: %f\n",(double)wid/returnwidcount());
  
/*As we will not sample from that distribution any longer, we free the allocated memory*/
  freesetup(t8);

  exit(0);
}

