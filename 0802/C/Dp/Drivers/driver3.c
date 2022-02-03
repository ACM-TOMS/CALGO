/* ALC2DLIB-Version 1.0 implemented by WH 22.3.99, all rights reserved */
/*please report problems or bugs to whoer@statistik.wu-wien.ac.at      */
/*                               or hormannw@boun.edu.tr               */

/*This file: mainnort.c */

#define PORTABLE 1 /*0..erf() is included in math-library, 1..erf() is not included*/

#define WID 100000
#define KLAS 30

#define START 17899234  /* starting value for uniform generator */


#define R 0.8
#define S1 1.
#define S2 1.




#include <math.h>
#include <stdio.h>
#include "alc2d.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

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
/*par[0]..R, par[1]..S1, par[2]..S2*/

double hn(double x[2],double par[]) /*log of bivariate normal density*/
{return
(-0.5/(1.-par[0]*par[0])*(x[0]*x[0]*(1./(par[1]*par[1]))-
    2.*par[0]*x[0]*x[1]*(1./(par[1]*par[2]))+x[1]*x[1]*(1./(par[2]*par[2]))));}

double hnx(double x[2],double par[]) /*Partielle Ableitung hn nach x*/
{return(-0.5/(1.-par[0]*par[0])*
	(2*x[0]*(1./(par[1]*par[1]))-2.*par[0]*x[1]*(1./(par[1]*par[2]))));}

double hny(double x[2],double par[]) /*Partielle Ableitung hn nach y*/
{ return(-0.5/(1.-par[0]*par[0])*
	 (-2.*par[0]*x[0]*(1./(par[1]*par[2]))+2*x[1]*(1./(par[2]*par[2]))));}

/*******************************************************************/
/*******************************************************************/
#if PORTABLE==0
double dfnormal(x) /*distribution function normal */
double x;
{
   if (x>=0) return((1+erf(x/M_SQRT2))/2);
   else return(erfc(-x/M_SQRT2)/2);
}
#else
/****************************************************************/
/* Verteilungsfunktion von N(0,1) programmiert nach             */
/* Kennedy/Gentle "Statistical Computing" S 90 - 92.            */
/****************************************************************/
double dfnormal(x)
double x;
{
  int  help;
  double xx,x2,x3,x4,x5,x6,x7,x8,zahler,nenner;

  static double
        ep0=242.66795523053175,ep1=21.979261618294152,
        ep2=6.9963834886191355,ep3= -3.5609843701815385e-2,
        eq0=215.0588758698612,eq1=91.164905404514901,
        eq2=15.082797630407787,/*eq3=1,*/
        zp0=300.4592610201616005,zq0=300.4592609569832933,
        zp1=451.9189537118729422,zq1=790.9509253278980272,
        zp2=339.3208167343436870,zq2=931.3540948506096211,
        zp3=152.9892850469404039,zq3=638.9802644656311665,
        zp4=43.16222722205673530,zq4=277.5854447439876434,
        zp5=7.211758250883093659,zq5=77.00015293522947295,
        zp6=0.5641955174789739711,zq6=12.78272731962942351,
        zp7= -1.368648573827167067e-7,/*zq7=1,*/
        dp0= -2.99610707703542174e-3,dq0=1.06209230528467918e-2,
        dp1= -4.94730910623250734e-2,dq1=1.91308926107829841e-1,
        dp2= -2.26956593539686930e-1,dq2=1.05167510706793207,
        dp3= -2.78661308609647788e-1,dq3=1.98733201817135256,
        dp4= -2.23192459734184686e-2/*,dq4=1*/;


  if(x<0)
    help=0;
  else
    help=1;
  x=fabs(x)/sqrt(2.);
  if(x<0.5)
    {
     x2=x*x;
     x4=x2*x2;
     x6=x4*x2;
     zahler=ep0+ep1*x2+ep2*x4+ep3*x6;
     nenner=eq0+eq1*x2+eq2*x4+x6;
     if(help)
       return(0.5*(1+x*zahler/nenner));
     else
       return(0.5*(1-x*zahler/nenner));
    }
  else if(x<4)
    {
     x2=x*x;
     x3=x2*x;
     x4=x3*x;
     x5=x4*x;
     x6=x5*x;
     x7=x6*x;
     zahler=zp0+zp1*x+zp2*x2+zp3*x3+zp4*x4+zp5*x5+zp6*x6+zp7*x7;
     nenner=zq0+zq1*x+zq2*x2+zq3*x3+zq4*x4+zq5*x5+zq6*x6+x7;

     if(help)
       return(0.5*(2-exp(-x2)*zahler/nenner));
     else
       return(0.5*exp(-x2)*zahler/nenner);
    }
  else if(x<50)
    {
     xx=x*x;
     x2=1/xx;
     x4=x2*x2;
     x6=x4*x2;
     x8=x6*x2;
     zahler=dp0+dp1*x2+dp2*x4+dp3*x6+dp4*x8;
     nenner=dq0+dq1*x2+dq2*x4+dq3*x6+x8;
     if(help)
       return(1-0.5*(exp(-xx)/x)*
             (1/sqrt(PI)+zahler/(xx*nenner)));
     else
       return(0.5*(exp(-xx)/x)*(1/sqrt(PI)+zahler/(xx*nenner)));
    }
  else if(help)
         return(1.);
       else
         return(0.);
}
#endif
/****************************************************************
 * cumulated distribution function - CHI^2 distribution         *
 * for n degrees of freedom,                                    *
 * Approximation from Wilson & Hilferty                         *
 * Kendall Stuart I, p 371                                      *
 ****************************************************************/
double dfchia(x,n)             /* approximation */
double x;
long n;
{
  double chix,y;

  y=2./(9.*n);
  chix=(exp(log(x/n)/3)-1+y)/sqrt(y);
  return(dfnormal(chix));
}
/****************************************/

/*Chi-2 test for uniform distribution*/ 

double chi2test(long b[]/*obsereved frequencies*/, 
                long l,/*number of classes*/
		long wid,/*sample size*/ 
		int printflag/*0..no output, 1.. little, 2..more*/)
{  long int i;
   double chi2=0.,erw,pval;

   erw=(double)wid/l;
   if (erw<5.) printf("Error chi2test: expected frequency smaller than 5\n");
   for (i=0;i<l;++i)
     {
       if (printflag==2) printf("%ld:%ld;",i,b[i]);
       chi2+= (b[i]-erw)*(b[i]-erw);
     }
   chi2/=erw;
   pval=1.-dfchia(chi2,l-1);
   if (printflag>=1)
     {
       printf(" Chi2-test: samplesize=%ld  number of classes= %ld \n",wid,l);
       printf("Chi2-value %f   Approximate P-Value:  %f\n",chi2,pval);
     }
   return(chi2);
}

/*******************************************************************/

long int bx[KLAS];
long int by[KLAS];
long int bw[KLAS];
long int bg[KLAS][KLAS];

#define HD 3.

int main()
{ long int j,ix,iy; 
  double x[2],z[2],sp[5][2],ineq[10][3],par[10],w;
/*array of the inequalities defining the domain */
/*inequalities are stored as 0 <= ineq[0][0] + ineq[0][1]*x + ineq[0][2]*y */  
  void *t1=NULL;


  ineq[0][0]=HD;
  ineq[0][1]=1.;
  ineq[0][2]=0.;
  ineq[1][0]=HD;
  ineq[1][1]=0.;
  ineq[1][2]=1.;
  ineq[2][0]=-HD;
  ineq[2][1]=1.;
  ineq[2][2]=0.;
  ineq[3][0]=-HD;
  ineq[3][1]=0.;
  ineq[3][2]=1.;

/*sets the seed for the uniform generator*/
   urn=START;

  sp[0][0]=0.1;
  sp[0][1]=0.2; 


  par[0]=R;par[1]=S1;par[2]=S2;
  printf("\nsetup and sample of size %d for normal-distr.: R=%f S1=%f S2=%f\n",WID,par[0],par[1],par[2]);

  t1=setup(1,20,sp,0,ineq,4,hn,hnx,hny,par,3);

  for(j=0;j<WID;j++)
  { 
    sample2d(x,t1);

/* The generated pair x is transformed into the standard normal pair z */
    z[0]=x[0]*(1./par[1]);
    z[1]=-x[0]*(par[0]/(par[1]*sqrt(1.-par[0]*par[0])))+x[1]*(1./(par[2]*sqrt(1.-par[0]*par[0])));

/* ix and iy are the number of the class, in which the generated pair falls */
    ix=(int)((dfnormal(z[0]))*KLAS);
    iy=(int)((dfnormal(z[1]))*KLAS);
/* the observed numbers are counted in:
           bx (marginal x) by (marginal y) bg (two dimensional) bw (angle) */
    ++bx[ix];
    ++by[iy];
    ++bg[ix][iy];
    w=atan2(z[1],z[0])+(PI);
    ++bw[(int)(w/(2.*PI)*KLAS)];
  }


  
  printf("Observed average acceptance probability: %f\n",(double)WID/returnwidcount());

  printf("volume below the hat approximately %f\n",
	 6.283185307*sqrt(S1*S1*S2*S2*(1.-R*R))/((double)WID/returnwidcount()));

  freesetup(t1);
  t1=NULL;
  printf("Marginal x:\n");
  chi2test(bx,KLAS,WID,1);
  printf("Marginal y:\n");
  chi2test(by,KLAS,WID,1);
  printf(":\n");
  chi2test(bw,KLAS,WID,1);
  printf("chi2 2 dim: \n");
  chi2test(bg[0],KLAS*KLAS,WID,1);

  exit(0);
}









