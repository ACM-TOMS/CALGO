/*------------------------------------------------------*/
/*  For use in C++Builder6 Under Windows */

/*
#include <vcl.h>
#pragma hdrstop
#pragma argsused
*/

/*------------------------------------------------------*/



/* New standard lib.*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <complex>
using namespace std;
typedef complex<double> complexn;

 /*Old standard */
/*
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <complex.h>
#define complexn  complex
*/



inline double acosh(double z)
       {return log(z+sqrt(z*z-1));}
inline double asinh(double z)
       {return log(z+sqrt(z*z+1));}

const double pi=3.141592653589793;


/*
  Copy Right Prof. Fayez A. Alhargan 2005.
  All Rights Reserved to the Author (Fayez A. Alhargan)
  No part or all of these subroutines may be used commercially
  without the written permission of the Author.
  Permission for academic use is given provided that due credit
  to the Author is given in any publications where these routines
  have been used for generating  some or all of the publication's results.
  The Author would be happy to help in any way he can
  in the use of these routines.

  Contact Address:
  CERI, KACST, P.O. Box 6086,
  Riyadh 11442,
  Saudi Arabia,
  Fax:966+1+4813764
  Email: alhargan@kacst.edu.sa

  Legal Matters:
  The program is here as is, no guarantees are given, stated or implied,
  you may use it at your own risk.

  Last updated 6 Nov 2005
*/

/*
  These routines are for testing:
  Mathieu Characteristic Numbers of general order
  and Mathieu Coefficients of general order

*/

      /*  Functions from the file:mcnf.cpp  */
  double MCNint(char typ,int n, double h);
	 /* Function used to obtain  c */
   double Compute_c(double v,double mu, double h,double &c2,double &acc);
	/* Function used to obtain real and imaginary  v */
   int  Compute_v(double c, double h, double &beta, double &mu,double &acc);
   int Approx_v(double c, double h, int &n, double &vre,double &vim);

    double CoefficientFr(double c, double v ,double h,
		     double pBem[],double nBem[],double pBom[],double nBom[],int CDIM, int MF);
    double CoefficientFr(complexn c, complexn v ,complexn h,
		     complexn pBem[],complexn nBem[],complexn pBom[],complexn nBom[],int CDIM, int MF);

                     
              /* Local functions for testing */
     double DetEc(FILE *val1,double c,double h);
     int Detminv(double a, double h, double Dsn, int &n, double &vre,double &vim);


int main()
{

   double h,q,t;
   int n,ne,i,vcmplx;
   double c,c2,cs,v,beta,mu,acc;
   double betaest,muest;
   double ae,bo,cl,cu;
   double pBem[200],nBem[200],pBom[200],nBom[200];
   complexn pBemx[200],nBemx[200],pBomx[200],nBomx[200];
   complexn cx,vx,hx,h2x,ux;

   double SB;
   int CDIM=50,MF;
   int m,mi,p;
   FILE *val1;
   double data[][4]={{1,0.5,1.0,2.537180071199},
		   {6,0.99,9.32355002165,49.7809026145324},
		   {10,0,0.1,100.0000505050675},
		   {10,0.00000001929,0.1,100.00005089086731},
		   {8,0.5,1000.0,-962.419932190115}}; /* Data{n,v,q,a} from shirts paper p.399-403*/
	       /* Coefficients Data from Shirts paper */
double TableV[][3]={{29,+1.0000000019e+01,+1.0000000000e-01},
{-18,+0.0000000000e+00,+2.5435707059e-33},
{-16,+5.6974599623e-30,-5.6975962214e-30},
{-14,-8.8880346066e-27,+8.8882446270e-27},
{-12,+8.5325029537e-24,-8.5327045730e-24},
{-10,-3.7542880298e-21,+3.7543767420e-21},
{ -8,-1.0457577300e-23,+1.0457823846e-23},
{ -6,-1.0445400278e-23,+1.0445196824e-23},
{ -4,-6.6743096113e-21,+6.6744734836e-21},
{ -2,-5.6066403340e-18,+5.6065506880e-18},
{ +0,-5.3822821095e-15,+5.3822848436e-15},
{ +2,-5.3822824482e-12,+5.3822819761e-12},
{ +4,-5.1669879094e-09,+5.1669880497e-09},
{ +6,-4.3402672010e-06,+4.3402672010e-06},
{ +8,-2.7777680404e-03,+2.7777680404e-03},
{+10,-9.9999355934e-01,+9.9999355934e-01},
{+12,+2.2727206201e-03,-2.2727206201e-03},
{+14,-2.3674185543e-06,+2.3674201351e-06},
{+16,+0.0000000000e+00,-1.5175779330e-09},
{+18,+0.0000000000e+00,+6.7749040131e-13},
{+20,+0.0000000000e+00,-2.2583019110e-16},
{+22,+0.0000000000e+00,+5.8809956481e-20},
{+24,+0.0000000000e+00,-1.2355034622e-23},
{+26,+0.0000000000e+00,+2.1449715287e-27},
{+28,+0.0000000000e+00,-3.1359235733e-31},
{+30,+0.0000000000e+00,+3.9199047633e-35},
{+32,+0.0000000000e+00,-4.2423214438e-39},
{+34,+0.0000000000e+00,+4.0173500726e-43},
{+36,+0.0000000000e+00,-3.3589885149e-47},
{+38,+0.0000000000e+00,+2.4992474988e-51}};


  val1=fopen("data.txt","wt") ;
      /* Generat table 1  */

   c=80.5;h=10.0;
   fprintf(val1,"\n Table 1");
   DetEc(val1,c,h);
   fprintf(val1,"  \\hline");
   printf("\n Finished Table 1");

	      /* Generat table 2  */
    ne=10;beta=0.5;mu=0;t=0;
    v=ne+beta;
    fprintf(val1,"\n Table 2");
     fprintf(val1,"  $\\nu=$%1.4lf",ne+beta);
   for(t=0;t<=4.0;t=t+0.25)
   {
    if(ne>0)
     {
      ae=MCNint('e',ne,t*ne)/(ne*ne);  /* lower limit is b_{n+1}/(n+1)^2 and upper limit is a_{n}/n^2 */
      bo=MCNint('o',ne+1,t*(ne+1.0))/((ne+1.0)*(ne+1.0));
      cl=bo+0.1*(1-beta)*sqrt(t)*(ae-bo);   /*  Equ. 28 */
      cu=ae-0.75*beta*(ae-bo);   /*  Equ. 29 */
      cl=cl*v*v;
      cu=cu*v*v;
    }
   else
    {
	/* for mu=0 and ne=0, lower limit is a_0 and upper limit is b_1 */
      ae=MCNint('e',ne,h);  /* equ(30) */
      bo=MCNint('o',ne+1,h); /* equ(31) */
      cu=bo;
      cl=ae;

    }
     h=t*v;
     c=Compute_c(ne+beta,mu,h,c2,acc);
     fprintf(val1,"\n");

     fprintf(val1,"  %1.2lf",t);
     fprintf(val1," & %1.6lf",cl);
     fprintf(val1," & %1.6lf",c);
     fprintf(val1," & %1.6lf",cu);
     fprintf(val1," & %1.4lf",c-cl);
     fprintf(val1," & %1.4lf",cu-c);
     fprintf(val1," \\\\ ");
   }
   fprintf(val1,"  \\hline");
  printf("\n Finished Table 2");

	    /* Generat table 3  */
    beta=0.5;mu=0;v=1.5;
    t=2.5;
    fprintf(val1,"\n Table 3");
     fprintf(val1,"  $t=$%1.4lf",t);
    for(v=1.5;v<101;v=v+10)
      {
	if(v==11.5) v=10.5;
	ne=int(v);
      if(ne>0)
      {
      ae=MCNint('e',ne,t*ne)/(ne*ne);  /* lower limit is b_{n+1}/(n+1)^2 and upper limit is a_{n}/n^2 */
      bo=MCNint('o',ne+1,t*(ne+1.0))/((ne+1.0)*(ne+1.0));
      cl=bo+0.1*(1-beta)*sqrt(t)*(ae-bo);   /*  Equ. 28 */
      cu=ae-0.75*beta*(ae-bo);   /*  Equ. 29 */
      cl=cl*v*v;
      cu=cu*v*v;
      }
     else
      {
	/* for mu=0 and ne=0, lower limit is a_0 and upper limit is b_1 */
      bo=MCNint('o',ne+1,h); /* equ(30) */
      ae=MCNint('e',ne,h);  /* equ(31) */
      cu=bo;
      cl=ae;
      }
     h=t*v;
     c=Compute_c(v,mu,h,c2,acc);
     fprintf(val1,"\n");
     fprintf(val1,"  %1.2lf",v);
     fprintf(val1," & %1.6lf",cl);
     fprintf(val1," & %1.6lf",c);
     fprintf(val1," & %1.6lf",cu);
     fprintf(val1," & %1.2lf",c-cl);
     fprintf(val1," & %1.2lf",cu-c);
     fprintf(val1," \\\\ ");
    }
   fprintf(val1,"  \\hline");
  printf("\n Finished Table 3");

	      /* Generat table 4  */
    ne=0;beta=0;mu=10.9;h=0.25;
    v=ne+beta;
    fprintf(val1,"\n Table 4");
     fprintf(val1,"  $\\nu=0+i$%1.4lf",mu);
   for(h=0.25;h<=101;h=h+10)
   {       /* this is for ne=0  only */
     ae=MCNint('e',ne,h);
     cu=ae;
     cl=ae-mu*mu*exp(-0.25*h*h/(mu*mu))-exp(-h/(4+mu*mu));
     c=Compute_c(v,mu,h,c2,acc);
     fprintf(val1,"\n");
     fprintf(val1,"  %1.2lf",h);
     fprintf(val1," & %1.6lf",cl);
     fprintf(val1," & %1.6lf",c);
     fprintf(val1," & %1.6lf",cu);
     fprintf(val1," \\\\ ");
   }
   fprintf(val1,"  \\hline");
  printf("\n Finished Table 4");

	      /* Generat table 5  */
    fprintf(val1,"\n Table 5");
    n=1;
    for(n=1;n<11;n+=1)
    {
     t=1; c=n*n*n+n*n+0.5;  h=n*3-2;
     fprintf(val1," \n %1.2lf",c);
     fprintf(val1," & %1.2lf",h);
	  /* Using the determinant method to estimate the order */
     vcmplx=Approx_v(c,h,ne,betaest,muest);
     if(vcmplx)  /*  i.e c lies in unstable region, v is complex  */
      {
       fprintf(val1," & %d+$i$%1.10lf",ne,muest);
       ne=Compute_v(c,h,beta,mu,acc);  /* Using bisection method */
       fprintf(val1," & %d+$i$%1.10lf",ne,mu);
       fprintf(val1," & $i$%1.2E",muest-mu);
       fprintf(val1," \\\\ ");
      }
     else   //  i.e c lies in stable region, v is real
      {
       mu=0;
       fprintf(val1," & %1.10f",ne+betaest);
       ne=Compute_v(c,h,beta,mu,acc);
       fprintf(val1," & %1.10f",ne+beta);
       fprintf(val1," & %1.2E",beta-betaest);
       fprintf(val1," \\\\ ");
      }
    }
   fprintf(val1,"  \\hline");
  printf("\n Finished Table 5");

	      /* Generat table 6  */
    fprintf(val1,"\n Table 6");
    n=1;
    for(n=1;n<11;n+=1)
    {
       t=3.5;  c=n*n+n+0.5;  h=t*n;  /* gives complex order  */
     fprintf(val1," \n %1.2lf",c);
     fprintf(val1," & %1.2lf",h);
	  /* Using the determinant method to estimate the order */
     vcmplx=Approx_v(c,h,ne,betaest,muest);
     if(vcmplx)  /*  i.e c lies in unstable region, v is complex  */
      {
       fprintf(val1," & %d+$i$%1.10lf",ne,muest);
       ne=Compute_v(c,h,beta,mu,acc);  /* Using bisection method */
       fprintf(val1," & %d+$i$%1.10lf",ne,mu);
       fprintf(val1," & $i$%1.2E",muest-mu);
       fprintf(val1," \\\\ ");
      }
     else   //  i.e c lies in stable region, v is real
      {
       mu=0;
       fprintf(val1," & %1.10f",ne+betaest);
       ne=Compute_v(c,h,beta,mu,acc);
       fprintf(val1," & %1.10f",ne+beta);
       fprintf(val1," & %1.2E",beta-betaest);
       fprintf(val1," \\\\ ");
      }
    }
   fprintf(val1,"  \\hline");
  printf("\n Finished Table 6");

      /* Generate c' table for  comparison between Shirts  and this method  */
     fprintf(val1,"\n Table 7");
  i=0;
  for(i=0;i<5;i++)
   {
     if(i==3) i++;
     ne=(int) data[i][0];beta=data[i][1];q=data[i][2];cs=data[i][3];
     h=2*sqrt(q);
     mu=0;
     c=Compute_c(ne+beta,mu,h,c2,acc);
     fprintf(val1,"\n");
     fprintf(val1,"  %1.4lf",ne+beta);
     fprintf(val1," & %1.6lf",h);
     fprintf(val1," & %1.13lf",cs);
     fprintf(val1," & %1.13lf",c);
     fprintf(val1," & %1.1E",c-cs);
//     fprintf(val1," & %1.1E",acc);
     fprintf(val1," \\\\ ");
   }
   fprintf(val1,"  \\hline");
  printf("\n Finished Table 7");

   /*  Generate Table 8 */
   fprintf(val1,"\n Table 8");
   MF=2;	 /*  McLachlan's  Normalization */
   i=3;
   ne=(int) data[i][0];beta=data[i][1];q=data[i][2];c=data[i][3];
   h=2.0*sqrt(q);
   mu=0;
   CoefficientFr(c,ne+beta,h,pBem,nBem,pBom,nBom,CDIM,MF);
   p=ne%2;
     fprintf(val1,"\n");
     fprintf(val1,"  $\\nu=$%1.11lf",ne+beta);
     fprintf(val1,"  $h=$%1.6lf ",h);
     fprintf(val1,"  $c=$%1.14lf",c);

   mi=(int) TableV[0][0];
   for(i=1;i<=mi;i++)
    {
     m=(int) TableV[i][0];
     SB=TableV[i][2];  /*  Shirts Coefficient */
     fprintf(val1,"\n");
     if(m<p)
	{
	  fprintf(val1," %3d ",m);
	  fprintf(val1," & %1.10E",SB);
	  fprintf(val1," & %1.10E",nBem[-m+p]);
	  fprintf(val1," & %1.2E",SB-nBem[-m+p]);
	  fprintf(val1," \\\\ ");

	}
      else
       {
	  fprintf(val1," %3d ",m);
	  fprintf(val1," & %1.10E",SB);
	  fprintf(val1," & %1.10E",pBem[m-p]);
	  fprintf(val1," & %1.2E",SB-pBem[m-p]);
	  fprintf(val1," \\\\ ");
       }
    }
	  fprintf(val1,"  \\hline");
  printf("\n Finished Table 8");
   /*  Generate Table 9 */
      fprintf(val1,"\n Table 9");
   MF=2;	 /*  McLachlan's  Normalization */
   h=2.0;
   mu=0.3;beta=0;
   ne=1;
   hx=complexn(h,0);
   h2x=0.25*hx*hx;
   vx=complexn(ne,mu);
   c=Compute_c(ne+beta,mu,h,c2,acc);
   cx=complexn(c,0);
   CoefficientFr(cx,vx,hx,pBemx,nBemx,pBomx,nBomx,CDIM,MF);
   p=ne%2;
     fprintf(val1,"\n");
     fprintf(val1,"  $\\nu=%1.4lf+i%1.4lf$",ne+beta,mu);
     fprintf(val1,"  $h=$%1.6lf ",h);
     fprintf(val1,"  $c=$%1.13lf",c);


   mi=15;
   m= -14;
   for(i=1;i<mi;i++)
    {
     m=m+2;
     fprintf(val1,"\n");
     if(m<0)
	{
	  ux=complexn(m+beta-p,mu);
	  acc=abs((nBemx[-m-2+p]+nBemx[-m+2+p])*h2x-(cx-ux*ux)*nBemx[-m+p]);
	  fprintf(val1," %3d ",m);
	  fprintf(val1," & %1.10E",real(nBemx[-m+p]));
	  SB=imag(nBemx[-m+p]);
	  if(SB<0) fprintf(val1,"  $-i$%1.10E",-SB); else fprintf(val1,"  $+i$%1.10E",SB);
	  fprintf(val1," & %1.2E",acc);
	  fprintf(val1," \\\\ ");

	}
    if(m==0)
	{
	  ux=complexn(m+beta+p,mu);
	  acc=abs((pBemx[2+p]+nBemx[2-p])*h2x-(cx-ux*ux)*pBemx[p]);
	  fprintf(val1," %3d ",m);
	  fprintf(val1," & %1.10E",real(nBemx[-m+p]));
	  SB=imag(nBemx[-m+p]);
	  if(SB<0) fprintf(val1,"  $-i$%1.10E",-SB); else fprintf(val1,"  $+i$%1.10E",SB);
	  fprintf(val1," & %1.2E",acc);
	  fprintf(val1," \\\\ ");
	}
    if(m>0)
       {
	  ux=complexn(m+beta+p,mu);
	  acc=abs((pBemx[m-2+p]+pBemx[m+2+p])*h2x-(cx-ux*ux)*pBemx[m+p]);
	  fprintf(val1," %3d ",m);
	  fprintf(val1," & %1.10E",real(pBemx[m-p]));
	  SB=imag(pBemx[m-p]);
	  if(SB<0) fprintf(val1,"  $-i$%1.10E",-SB); else fprintf(val1,"  $+i$%1.10E",SB);
	  fprintf(val1," & %1.2E",acc);
	  fprintf(val1," \\\\ ");
       }
    }
     ux=complexn(beta+p,mu);
     acc=abs((pBemx[2+p]+nBemx[2-p])*h2x-(cx-ux*ux)*pBemx[p]);

    fprintf(val1,"  \\hline");
  printf("\n Finished Table 9\n");

/*=========================================================*/
//   cin>>mi;
   fclose(val1); exit(0);

}

double DetEc(FILE *val1,double c,double h)
{
      /*
	Here v=0.
	This function computes  E(c,h).
	E(c,h)=Dlt(c,0,h)*sin^2(0.5*pi*c^0.5) given by McLachlan p.69.
	Also see companion paper for Mathematical derivation,
	note singularities when c nears 4*n*n;  n=0,1,2...
	these are taken care of.
      */
  double E0,E1,E2;
  double h2,acc=1,csr;
  double sn,sn2;
  double b0,b1,b2,bn,d0,d1,d2,dn;
  double bm1,Em1,Em2;
  int RL,m,n=-1;
  int nr,r; double beta,vim;

  h2=h*h*0.25;
     fprintf(val1," \n & %1.5lf",c);
     fprintf(val1," & %1.5lf",h);

  if(c>=0){
	     csr=sqrt(c);
	     sn=sin(0.5*pi*csr);
	     sn2=sn*sn;
	     n=int(0.5*csr);
	     if(fabs(sn)<1.0e-4)
		   {
		     if(n==0) {bn=-h2*pi*pi*0.25;dn=sn2;}
		       else {bn=-h2*pi*cos(0.5*pi*csr)/(4.0*csr);dn=sn;}
		    }
		 else
		    {
		      if(n==0) {bn=h2/(4.0*n*n-c)*sn2;dn=sn2;}
			 else {bn=h2/(4.0*n*n-c)*sn;dn=sn;}
		    }

	  }
	 else
	  {            /* c<0,  no singularity   */
	     /* sin term has been taken outside and mutliplied at the end */
	     csr=sqrt(-c);
	     sn=sinh(0.5*pi*csr);
	     sn2=-sn*sn;
	     bn=h2/(4.0*n*n-c);
	     dn=1;
	  }

  RL=8000;
  if(RL<3*fabs(c)) RL=(int)(5.0*fabs(c)+20);

  d0=1;d1=1;d2=1;
  if(n==0)  {b0=bn;d0=dn;} else  {b0=-h2/c;d0=1;}
  if(n==1)  {b1=bn;d1=dn;} else  {b1=h2/(4.0-c);d1=1;}
  if(n==2)  {b2=bn;d2=dn;} else  {b2=h2/(4.0*2*2-c);d2=1;}

  Em2=0;
  Em1=0;
  E0=d0;
  E1=d0*d1*d1-2*b0*b1*d1;

   Detminv(c,h,E1,nr,beta,vim);
   fprintf(val1," \n %d & %1.8lf",1,E1);
   fprintf(val1,"   & %1.8lf",beta);
   if(vim!=0)  fprintf(val1,"+$i$%1.8lf",vim);
   fprintf(val1," \\\\ ");

  E2=d0*pow(d2*d1-b2*b1,2)-2*d2*b1*b0*(d2*d1-b2*b1);
   Detminv(c,h,E2,nr,beta,vim);
   fprintf(val1," \n %d & %1.8lf",2,E2);
   fprintf(val1,"   & %1.8lf",beta);
   if(vim!=0)  fprintf(val1,"+$i$%1.8lf",vim);
   fprintf(val1," \\\\ ");

  m=2;r=2;
  while(m<n+4)
    {
	  m++;
	  Em2=Em1;
	  Em1=E0;
	  E0=E1;
	  E1=E2;
	  bm1=b0;
	  b0=b1;d0=d1;
	  b1=b2;d1=d2;
	  if(m==n)  {b2=bn;d2=dn;}  else  {b2=h2/(4.0*m*m-c);d2=1;}

	  if(fabs(d0)<0.1)
	      {
		  E2=(d2*d1-b2*b1)*d2*E1-d1*b2*b1*(d2*d1-b2*b1)*E0
			+d2*b0*b0*b1*b1*b1*b2*Em1;
		  E2=E2/d1;
	      }
	    else
	     {
	       E2=d2*d2*d0*E1-d0*b2*b1*(2*d2*d1-b2*b1)*E0+d2*b2*b1*b1*b0*E0
		+d2*d0*d0*b2*b1*b1*b0*Em1-d2*b2*b1*b1*b0*b0*b0*bm1*bm1*Em2;
	       E2=E2/d0;
	      }
	r++;
//	if((r+4)%7==0)
         if((log(1.0*r)/log(2.0)- int(log(1.0*r)/log(2.0)))==0)
	    {
	     Detminv(c,h,E2,nr,beta,vim);
	     fprintf(val1," \n %d & %1.8lf",r,E2);
	     fprintf(val1,"   & %1.8lf",beta);
 	     if(vim!=0)  fprintf(val1,"+$i$%1.8lf",vim);
	     fprintf(val1," \\\\ ");
	    }

	  acc=fabs(E1-E2);
     }

  while(acc>1e-18 && m<RL+2)
    {
	  m++;
	  Em2=Em1;
	  Em1=E0;
	  E0=E1;
	  E1=E2;
	  bm1=b0;
	  b0=b1;
	  b1=b2;
	  b2=h2/(4.0*m*m-c);
	  E2=E1-b2*b1*(E1+(1-b2*b1)*E0-b0*b0*b1*b1*Em1);
	 r++;
//	if((r+4)%7==0)
         if((log(1.0*r)/log(2.0)- int(log(1.0*r)/log(2.0)))==0)
	    {
	     Detminv(c,h,E2,nr,beta,vim);
	     fprintf(val1," \n %d & %1.8lf",r,E2);
	     fprintf(val1,"   & %1.8lf",beta);
	    if(vim!=0)  fprintf(val1,"+$i$%1.8lf",vim);
	    fprintf(val1," \\\\ ");
	    }

	  acc=fabs(E1-E2);
     }

   if(c<0) E2=E2*sn2;


   return E2;
}
int Detminv(double c, double h, double Dsn, int &n, double &vre,double &vim)
{
  /*
    This is a test routine only
    Given c and h find approximate value of v.
    1. Locate the position of c on the h-c plan to get the value of n
    2. Evaluate the Determinant
    3. Find vre and vim using the standard  equation
    4. return 0 if v is real and  1 if v is complex
  */

 double v,t,ae,bo;
 int ne,vcmplx;



   n=0;vre=0;vim=0;

   ne=0;vcmplx=0;
   ae=MCNint('e',ne,h);
   bo=MCNint('o',ne+1,h);
   while(!(c>ae && c<bo))
   {
     ne++;
     ae=MCNint('e',ne,h);
     bo=MCNint('o',ne+1,h);
     if(c<ae) {vcmplx=1;break;}
   }

   if(vcmplx==1)  /*  i.e a lies in unstable region */
   {       /* Find the integer part of the order v */
	ne=0;
	ae=MCNint('e',ne,h);
	if(ne>0) bo=MCNint('o',ne,h); else bo=-1e200;
	while(!(c>bo && c<ae))
	{
	   ne++;
	   ae=MCNint('e',ne,h);
	   bo=MCNint('o',ne,h);
	}
   }


   if(Dsn<0)
	 {
	   t=sqrt(-Dsn);
	   vim=2*asinh(t)/pi;   /* imaginar, unstable region */
	   vre=0;
	   n=ne;
	   vcmplx=1;
	 }
      else
	 {
	   t=sqrt(Dsn);
	   if (t>1)
	      {   /* imaginar, unstable region */
		if(ne%2==0)
		  vim=2*asinh(t)/pi;
		else
		  vim=2*acosh(t)/pi;
		vre=0;n=ne;
		vcmplx=1;
	      }
	     else
	      {
		v=2*asin(t)/pi;   /* stable region */
		if(ne%2==0)  vre=v; else vre=1-v;
		n=ne;
		vim=0;
		vcmplx=0;
	      }
	 }

 return vcmplx;

}
