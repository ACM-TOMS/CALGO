
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





const double pi=3.141592653589793;
inline double acosh(double z)
       {return log(z+sqrt(z*z-1));}
inline double asinh(double z)
       {return log(z+sqrt(z*z+1));}

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
  These routines compute:
  Mathieu Characteristic Numbers of general order
  and Mathieu Coefficients of general order

*/

/*--------------------------------------------------------------*/
/*                  Functions required from MCNR.CPP            */
/*--------------------------------------------------------------*/

     double Estimatmcn(char typ,int n, double h);
     double MCNRoot(double a, char typ, int n, double h, double accr, double &acc);


/*--------------------------------------------------------------*/
/*                     MCNF header file                         */
/*--------------------------------------------------------------*/
      /*
	 declears a pointer to a function
	 this minimize rewriting the bisection routine for each
	 function of f(c), f(v) or f(mu)
      */
  typedef double (*ptf)(double c,double v,double mu, double h);

       /*  complex function used in root finding routines */
   complexn CFmcnComp(double c,complexn v, double h);
	 /* real function used in root finding routines */
   double CFmcnReal(double c,double v, double h);

	  /* function for root finding of MCN c */
   double CFmcnfc(double c, double v,double mu, double h);
	  /* function for root finding of the real order v */
   double CFmcnfv(double v, double c,double mu, double h);
	 /* function for root finding of  complex order v+i mu */
   double CFmcnfmu(double mu,double c,double v, double h);

	   /* Search a function f: forward */
   void Searchf(ptf f ,double y, double z,double h, double x1,double sp,double &xs,double &xl);
	 /* Search a function f: backward */
   void Searchbk(ptf f ,double v,double mu, double h, double sp,double &xs,double &xl);
	  /* Use bisection to find a root of the function f */
   double Bisection(ptf f,double y,double z, double h, double xl,double xu, double &acc);

	  /* Compute Integer order  MCN  */
   double MCNint(char typ,int n, double h);

	 /* Function used to obtain  c */
   double Compute_c(double v,double mu, double h,double &c2,double &acc);
	/* Function used to obtain real and imaginary  v */
   int  Compute_v(double c, double h, double &beta, double &mu,double &acc);

       /*  Routine for computing the determinant E(c,h)*/
   double DetE(double c,double h);
       /* Routine used to obtain approximate value of v=n+beta+i mu */
   int Approx_v(double c, double h, int &n, double &beta,double &mu);

       /* Routine used to compute the Mathieu coefficients
		    of fractional order  */
   double CoefficientFr(double c, double v ,double h,
	     double pBe[],double nBe[],double pBo[],double nBo[],int CDIM, int MF);
   double CoefficientFr(complexn c, complexn v ,complexn h,
		     complexn pBem[],complexn nBem[],complexn pBom[],complexn nBom[],int CDIM, int MF);


/*===================================================================*/

/*===================================================================*/
double CFmcnfc(double c, double v,double mu, double h)
{
	/* Function used to evluate the continued fractions
	     as a function of c */
	 if (mu==0)
	    return CFmcnReal(c,v,h);
	 else
	    return real(CFmcnComp(c,complexn(v,mu),h));
}

double CFmcnfv(double v, double c,double mu, double h)
{
	/* Function used to evluate the continued fractions
	     as a function of v */

	 if (mu==0)
	    return CFmcnReal(c,v,h);
	 else
	    return real(CFmcnComp(c,complexn(v,mu),h));
}

double CFmcnfmu(double mu,double c,double v, double h)
{
	/* Function used to evluate the continued fractions
	     as a function of mu (imaginary part of the order) */

	 if (mu==0)
	    return CFmcnReal(c,v,h);
	 else
	   return real(CFmcnComp(c,complexn(v,mu),h));
}
/*===================================================================*/
/*

   Compute integer order MCN
   Input: typ:  'e':even  'o':odd
	    n:  order integer
	    h:  parameter

   Output: return integer order MCN

   Note the body of the routines used are in the file MCNR.CPP
*/
double MCNint(char typ,int n, double h)
{
    /* Evaluate integer order mcn available in mcnr.cpp */
    double R,acc;

    R=Estimatmcn(typ,n,h);  /* Estimate MCN of order n */
    acc=1e-8;
    R=MCNRoot(R,typ,n,h,1e-5,acc); /* Improve the Accuracy of MCN */

   return R;
}

/*====================================================================*/
/*---------------------------------------------------------------*/
double Compute_c(double vr,double mu, double h,double &c2,double &acc)
{
   /*
     vr=n+beta real part, mu imaginary part of the order
     Given v=vr+i mu and h, find c.
     1. Find lower and upper limits for c
     2. Use bisection method to improve accuracy

   */
    int n;
   double beta,t;
   double c,cl,cu,ae,bo;
   double cl2,cu2;
   double cl3,cu3;
   double sp;
   int ne;
   ptf f;

    n=(int) floor(vr+1.0e-9);  /* v very close to integers are approximately integers */
    beta=vr-n;

   f=CFmcnfc;  /* Pointing to CFmcnfc: function of c */

    vr=n+beta;
    ne=(int)floor(vr+1.0e-9);  /* v very close to integers are approximately integers */
     if(ne!=0) t=h/vr; else t=h;
  if(mu==0)  /* i.e. real order  */
  {
    if(fabs(vr-ne)<1e-9)
      {   /* it is an integer order MCN */
	   if(beta<0.5)
	      c=MCNint('e',ne,h);
	   else
	      c=MCNint('o',ne,h);
      }
   else
    {
    if(ne>0)
     {
      ae=MCNint('e',ne,t*ne)/(ne*ne);
      bo=MCNint('o',ne+1,t*(ne+1.0))/((ne+1.0)*(ne+1.0));

	  /* lower limit is b_{n+1}/(n+1)^2 and
	      upper limit is a_{n}/n^2 */
      cl=bo+0.1*(1-beta)*sqrt(t)*(ae-bo);   /*  Equ. 28 */
      cl=cl*vr*vr;
      cu=ae-(0.75*beta)*(ae-bo);   /* Equ. 29 */
      cu=cu*vr*vr;
      if(t>0 && t<=0.5) {cu=ae*vr*vr;cl=bo*vr*vr;}
    }
   else
    {
      /* for ne=0, lower limit is a_0 and upper limit is b_{1} */
      ae=MCNint('e',ne,h);
      bo=MCNint('o',ne+1,h);
      cu=bo;
      cl=ae;
    }

     sp=fabs(cu-cl)/40.0;
     if(sp==0) {sp=1.0e-9;cl=cl-sp;cu=cu+sp;}
    if(h==0)
    {
      c=vr*vr;
    }
   else
    {
      /*takes cu as upper and searches for improved lower&upper limits */

     Searchbk(f,vr,mu,h,sp,cl,cu);
     c=Bisection(f,vr,mu,h,cl,cu,acc);   /* bisection method */
    }
   }
   c2=c;
  }
 else
  {    /* complex order, two valued function  */
   if(ne>0)
    {
     ae=MCNint('e',ne,h);
     bo=MCNint('o',ne,h);
	  /* lower limit is b_{n} and upper limit is a_{n} */
     cl=bo;
     cu=ae;
     sp=fabs(cu-cl)/80.0;
     if(sp==0) {sp=1.0e-4;cl=cl-4.0*sp;cu=cu+4.0*sp;}
    }
   else
    {               /*  ne=0  */
     ae=MCNint('e',ne,h);
     cu=ae;
     cl=ae-mu*mu*exp(-0.25*h*h/(mu*mu))-exp(-h/(4+mu*mu));
     if(fabs(ae-cl)<0.1) cu=cu+0.01;
     sp=fabs(cu-cl)/40.0;
     if(sp<1e-2)
	{sp=1.0e-1;cl=cu-sp;cu=cu+sp;}
    }
    if(h==0)
    {
      c=vr*vr;
      c2=c;
    }
   else
    {
     cl2=cl;cu2=cu;
     cl3=cl;cu3=cu;
     Searchbk(f,vr,mu,h,sp,cl,cu);
     c=Bisection(f,vr,mu,h,cl,cu,acc);   /* bisection method  */
      if(c>cu3 || c<cl3) {cout<<"\n Error: Compute_c(),  obtained root outside the limits, vr="<<vr;}
   if(ne>0)
     {
      Searchf(f,vr,mu,h,cl2,sp,cl2,cu2);
      c2=Bisection(f,vr,mu,h,cl2,cu2,acc);    /* bisection method */
      if(c2>cu3 || c2<cl3) {cout<<"\n Error2: Compute_c(),  obtained root outside the limits, vr="<<vr;}
     }
    else
     c2=c;

    }

  }
   return c;

}

/*----------------------------------------------------*/
int  Compute_v(double c, double h, double &beta, double &mu,double &acc)
{
  /*
    Given c and h find the order v.
    1. Find approximate beta & mu  using determinant method
    2. Refine the the roots using CF method
    3. return n integer part
    4. return beta:fractional part,  mu:imaginary part
       and acc:accuracy of the computation

  */

 double v;
 double vl,sp,vu;
 int ne,vcmplx;


	    /* Find an approximate value of v */
   vcmplx=Approx_v(c,h,ne,beta,mu);v=ne+beta;

   if(vcmplx==1)  /*  i.e c lies in unstable region, v is complex */
   {
	beta=0;
	mu=mu-0.1;
	if(mu<0) mu=1e-9;
	   /* Search forward and  Bisect the function CFmcnfmu */
	Searchf(CFmcnfmu,c,ne,h,mu,0.1,vl,vu);
	mu=Bisection(CFmcnfmu,c,ne,h,vl,vu,acc);
   }
  else  /*  i.e c lies in stable region, v is real  */
   {
     mu=0;
     vl=(1-0.02/(fabs(c)+1))*v;
     vu=(1+0.02/(fabs(c)+1))*v;sp=fabs(vu-vl)/10.0;
	  /* Make sure vl,vu are within the reasonable limits */
     if(vu>(ne+1.0)) vu=ne+1.0;
     if(vl<ne)  vl=ne;
	       /* Search and Bisect the function CFmcnfv */
     Searchf(CFmcnfv,c,mu,h,vl,sp,vl,vu);
     v=Bisection(CFmcnfv,c,mu,h,vl,vu,acc);
     beta=v-ne;mu=0;
   }

 return ne;

}
/*---------------------------------------------------------*/
/*		Determinant Evaluation E(a,h)	           */
/*---------------------------------------------------------*/
double DetE(double c,double h)
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

  h2=h*h*0.25;

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

  RL=100;
  if(RL<3*fabs(c)) RL=(int)(5.0*fabs(c)+20);

  d0=1;d1=1;d2=1;
  if(n==0)  {b0=bn;d0=dn;} else  {b0=-h2/c;d0=1;}
  if(n==1)  {b1=bn;d1=dn;} else  {b1=h2/(4.0-c);d1=1;}
  if(n==2)  {b2=bn;d2=dn;} else  {b2=h2/(4.0*2*2-c);d2=1;}

  Em2=0;
  Em1=0;
  E0=d0;
  E1=d0*d1*d1-2*b0*b1*d1;
  E2=d0*pow(d2*d1-b2*b1,2)-2*d2*b1*b0*(d2*d1-b2*b1);
  m=2;
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

	  acc=fabs(E1-E2);
     }
   if(c<0) E2=E2*sn2;


   return E2;
}
/*---------------------------------------------------------*/
int Approx_v(double c, double h, int &n, double &beta,double &mu)
{
  /*
    Given c and h find approximate value of the order v=n+beta+ i mu.
    1. Locate the position of c on the h-c plan to get the value of n
    2. Evaluate the Determinant
    3. Find beta and mu using the standard  equation
    4. return 0 if v is real and  1 if v is complex
  */

 double v,Es,E,ae,bo;
 int ne,vcmplx;

   n=0;beta=0;mu=0;

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

   if(vcmplx==1)  /*  i.e c lies in unstable region */
   {
	 /* Find the integer part of the order v */
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

   n=ne;
   E=DetE(c,h);
   if(E<0)
	 {
	   Es=sqrt(-E);
	   mu=2*asinh(Es)/pi;   /* imaginary, unstable region */
	   beta=0;
	   vcmplx=1;
	 }
      else
	 {
	   Es=sqrt(E);
	   if (Es>1)
	      {   /* imaginary, unstable region */
		if(ne%2==0)
		  mu=2*asinh(Es)/pi;
		else
		  mu=2*acosh(Es)/pi;

		beta=0;
		vcmplx=1;
	      }
	     else
	      {
		v=2*asin(Es)/pi;   /* stable region */
		if(ne%2==0)  beta=v; else beta=1-v;
		mu=0;
		vcmplx=0;
	      }
	 }

 return vcmplx;

}

/*================================================================*/
/*    Evaluating MCNs using Continued Fractions			  */
/*================================================================*/
double CFmcnReal(double c, double v, double h)
{
    /* Continued Fractions function  see the companion paper */
   /* also see McLachlan p.107    */
   /* Singularity dealt with */
   /* Input: c,v,h   Output: fr */

   double fn,fp;
   double h2,R;
   double u,bt;
   double n;
   double Vm,Vm1,Ti,cr;
   long vi,m,K;
   int p;

    vi=int(v);
    bt=v-vi;
    p=vi%2;
    if(c>=0) cr=0.5*(sqrt(c)-bt-p); else cr=0;
    n=int(cr);
   h2=0.25*h*h;

   K=(int)fabs(c);
     /* computing for negative m */
  Vm1=Vm=0;
  for(m=40+K;m>0;m--)  /* backward */
    {
      u=(-2.0*m+bt+p);
      Ti=1.0/(u*u);
      Vm1=-Ti*h2/(1-c*Ti+Ti*h2*Vm);
      Vm=Vm1;
    }

   u=(bt+p);
   fn=(c-u*u)/h2-Vm1;

     /* computing for positive m */

   Vm1=fn;
  for(m=1;m<=n;m++)  /* forward */
   {
    u=(2.0*m+bt+p);
    Vm=(c-u*u)/h2-1.0/Vm1;  /* Vm1=V_{m-1}  and Vm=V_{m} */
    Vm1=Vm;
   }
   fn=Vm1;

  Vm=0;
  for(m=40+K;m>n;m--)  /* backward  */
    {
      u=(2.0*m+bt+p);
      Ti=1.0/(u*u);
      Vm1=-Ti*h2/(1-c*Ti+Ti*h2*Vm);
      Vm=Vm1;
    }
   fp=Vm;

   R=(fp-fn);  /* equation (13) in the companion paper */

   return R;
}
/*----------------------------------------------------------*/
complexn CFmcnComp(double c, complexn v, double h)
{
   /*
     same as the function CFmcnReal(c,v,h) above with order v being complex.
   */

   complexn fn,fp,R;
   complexn u,y,bt;
   complexn V0,Vm,Vm1,Ti;
   double h2;
   double n;
   double cr;
   long vi,m,K;
   int p;

    vi=int(real(v));
    bt=v-complexn(vi,0);
    p=vi%2;
    if(c>=0) cr=0.5*(sqrt(c)-real(bt)-p); else cr=0;
    if(imag(v)==0) n=int(cr); else  n=0;
    h2=0.25*h*h;
    K=(int)fabs(c);
     /* computing for negative m */
  Vm1=Vm=0;
  for(m=40+K;m>0;m--)  /* backward  */
    {
      u=complexn(-2.0*m+p,0)+bt;
      Ti=1.0/(u*u);
      Vm1=-Ti*h2/(complexn(1,0)-c*Ti+Ti*h2*Vm);
      Vm=Vm1;
    }
   u=(bt+complexn(p,0));
   fn=(c-u*u)/h2-Vm1;

     /* computing for positive m */
   Vm1=Vm=fn;
  for(m=1;m<=n;m++)  /* forward  */
   {
    u=complexn(2.0*m+p,0)+bt;
    Vm=c/h2-u*u/h2-1.0/Vm1;  /* Vm1=V_{m-1}  and Vm=V_{m} */
    Vm1=Vm;
   }
   fn=Vm1;

  Vm=0;
  for(m=40+K;m>(int)n;m--)  /* backward  */
    {
      u=complexn(2.0*m+p,0)+bt;
      Ti=1.0/(u*u);
      Vm1=-Ti*h2/(complexn(1,0)-c*Ti+Ti*h2*Vm);
      Vm=Vm1;
    }
   fp=Vm;

   R=(fp-fn);

   return R;
}
/*-----------------------------------------------------------------------*/
/*--------------- Find the lower and upper limits --------------------*/
void Searchbk(ptf f,double y,double z, double h, double sp,double &xl,double &xu)
{
  /*
    This is a searching algorithm which given an initial value(xi)
    and step (sp) will locate the root and gives xl and xu.
    This routine works backwards, i.e. xi>root.
    xi=xu=starting point, sp=step size, xl=lower value, xu=upper value.
  */
  double Fx1,Fx2,x1,x2;
  int m;

  m=0;x2=xu;x1=xu-sp;
  Fx1=f(x1,y,z,h);
  Fx2=f(x2,y,z,h);
  while (Fx2*Fx1>0  && m<150)
  {
   x2=x1;
   x1-=sp;
   Fx2=Fx1;
   Fx1=f(x1,y,z,h);
   m++;
  }
  if(Fx2*Fx1>0)
	{ cout << "\n Error: Searchbk() Finding Limits has failed y="
	      <<y<<" h="<<h;
	}
     else
      {
	xl=x1; xu=x2;
      }
}
/*--------------- Find the lower and upper limits --------------------*/
void Searchf(ptf f,double y, double z, double h, double xi,double sp,double &xl,double &xu)
{
  /*
    This a searching algorithm which given an initial value(xi)
    and a step (sp) will locate the root and obtains xs and xu.
    This routine works forwards, i.e. xi<root.
    xi=starting point, sp=step size, xl=lower value, xu=upper value.
  */
  double Fx1,Fx2,x1,x2;
  int m;

  m=0;x1=xi;x2=xi+sp;
  Fx1=f(x1,y,z,h);
  Fx2=f(x2,y,z,h);
  while (Fx1*Fx2>0  && m<150)
  {
   x1=x2;
   x2+=sp;
   Fx1=Fx2;
   Fx2=f(x2,y,z,h);
   m++;
  }
  if(Fx2*Fx1>0)
     cout << "\n Error: Searchf() Finding Limits has failed y="<<y;
   else
    {  xl=x1; xu=x2;}
}
/*------------ Bisection Method ------------------------------------*/
double Bisection(ptf f,double y,double z, double h, double xl,double xu,double &acc)
{
  /*
    This is a general root finding algorithm which uses
    bisection method to refine the root,
    given  lower (xl) and upper(xu) limits to get more accurate root.
    Ridders' Method is used, see numerical receipes in C p.358
  */
 double  x1,x2,x3,Fx1,Fx2,Fx3,m,tst,tole=1e-15;
 double  x4,Fx4,sq;
 int sign1,sign2,sign3,sign4;

  x1=xl;x2=xu;x3=x1;tst=1;m=0;
  x4=0.5*(x1+x2);
  Fx1=f(x1,y,z,h);
  Fx2=f(x2,y,z,h);
  if(Fx1<0) sign1=1; else sign1=-1;
  if(Fx2<0) sign2=1; else sign2=-1;

  acc=Fx1;
  while (fabs(acc)>=tole && m<100 && tst>tole)
  {
    x3=0.5*(x2+x1);
    Fx3=f(x3,y,z,h);
    if(Fx3<0) sign3=1; else sign3=-1;

    sq=sqrt(fabs(Fx3*Fx3-Fx1*Fx2));
    if(Fx3==0 || sq==0)  return x3;
    if((Fx1-Fx2)>=0)
	x4=x3+(x3-x1)*Fx3/sq;
      else
	x4=x3-(x3-x1)*Fx3/sq;

    Fx4=f(x4,y,z,h);
    if(Fx4<0) sign4=1; else sign4=-1;

    if(sign3!=sign4)
     {
      x1=x3;Fx1=Fx3;sign1=sign3;
      x2=x4;Fx2=Fx4;sign2=sign4;
     }
    else if(sign1!=sign4)
     {
      x2=x4;Fx2=Fx4;sign2=sign4;
     }
    else if(sign2!=sign4)
     {
      x1=x3;Fx1=Fx3;sign1=sign3;
     }
    else {cout<<"Error: function Bisection() is unable to get the root";exit(0);}

    tst=fabs(x1-x2);
    acc=Fx4;
    m++;
   }
  return x4;
}
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
double CoefficientFr(double c, double v ,double h,
		     double pBem[],double nBem[],double pBom[],double nBom[],int CDIM, int MF)
{
    /*
      You can chose several normalization
      0: No normalization
      1: Morse and Feshbach
      2: McLachlan (p.112)

      For coefficient computation
      inputs:
      c=MCN
      v=fractional order
      h=pamarater
      CDIM= Coefficients Array dimension

      outputs:
      pBem[]=Array of Coefficients of positive index for even functions
      nBem[]=Array of Coefficients of negative index for even functions
      pBom[]=Array of Coefficients of positive index for odd  functions
      nBom[]=Array of Coefficients of negative index for odd  functions
      local:
      pW= positive indexed array
      nW= negative indexed array
    */
    double *nW,*pW;
    double h2,hi2,bet,cr,acc;
    double nC1;
    int Ml,p,i,vi,r,m;
    int index;
    double Norme,Normo;
    double Vm,Vm1,Ti,u;

    vi=int(v);
    bet=v-vi;
    p=vi%2;


    if(c>=0) cr=0.5*(sqrt(c)-bet-p); else cr=0;
    r=int(cr);  /* Center the computation around the singularity */
    if(r<0) r=0;
    Ml=5+r+CDIM;  // Large m

   h2=0.25*h*h;
   hi2=1.0/h2;
    nW=new double [2*CDIM+6];
    pW=new double [2*CDIM+6];
    if(nW==NULL)
      {
       cout<<"\n CoefficientFr: out of memory..."; exit(0);
      }
    if(pW==NULL)
      {
       cout<<"\n CoefficientFr: out of memory..."; exit(0);
      }

    for(i=0;i<2*CDIM+2;i++)
     {
       nW[i]=0;
       pW[i]=0;
     }

     /* computing for negative m */
  Vm1=Vm=0;
  for(m=Ml;m>1;m--)  // backward  down to m=2
    {
      u=(-2.0*m+bet+p);
      Ti=1.0/(u*u);
      Vm1=-Ti*h2/(1-c*Ti+Ti*h2*Vm);
      Vm=Vm1;
      if(m<=CDIM+1) nW[2*m-p-2]=Vm1;
    }
   m=1;              // backward    m=1
   u=(-2.0*m+bet+p);
   Ti=1.0/(u*u);
   Vm1=-Ti*h2/(1-c*Ti+Ti*h2*Vm);
   nC1=Vm1;

     /* computing for positive m */
   Vm=(c-(bet+p)*(bet+p))/h2-Vm1;  // forward  m=0;
   Vm1=Vm;
   pW[p]=Vm;
  for(m=1;m<=r;m++)  // forward up to m=r
   {
    Ti=(2.0*m+bet+p)*(2.0*m+bet+p);
    Vm=c*hi2-Ti*hi2-1.0/Vm1;  // Vm1=V_{m-1}  and Vm=V_{m}
    Vm1=Vm;
    if(m<=CDIM) pW[2*m+p]=Vm;
   }

  Vm=0;
  for(m=Ml;m>r;m--)  // backward  down to m=r+1
    {
      u=(2.0*m+bet+p);
      Ti=1.0/(u*u);
      Vm1=-Ti*h2/(1-c*Ti+Ti*h2*Vm);
      Vm=Vm1;
      if(m<=CDIM) pW[2*m+p-2]=Vm;
    }

      pBem[p]=1;
      for(m=0;m<=CDIM;m++)
       {pBem[2*m+2+p]=pW[2*m+p]*pBem[2*m+p];}

      nBem[0]=0;
      nBem[2-p]=pBem[p]*nC1;
      for(m=1;m<=CDIM;m++)
       {nBem[2*m+2-p]=nW[2*m-p]*nBem[2*m-p];}

	 /* accuracy check see McLachlan p.113 eq(19)  */
	 /*     c1=pBm[p];c3=pBem[2+p]; cm1=nBem[2-p];  */
     u=bet+p;
     acc=((pBem[2+p]+nBem[2-p])*h2-(c-u*u)*pBem[p]);
     acc=fabs(acc);
     if(acc>1e-14) {cout<<"\n Warning: CoefficientFr() is not accurate acc="<<acc;}

   /*----------------------------------------------------------------*/

       /* Normalization convensions */
    if(MF==1)
    {      /*  Morse and Feshbach Normalization convension
	       in MF the even is different from the odd normalization
	   */
     index=p;  // i=0
     if(p==0)
     {
	Norme=pBem[index];
	Normo=(p+bet)*(pBem[index]);
     }
   else
    {
	 Norme=pBem[index]+nBem[index];  // add C[-1] and C[+1]
	 Normo=(p+bet)*pBem[index]+(p+bet)*nBem[index];
    }
     for(i=1;i<=CDIM;i++)
	{
	 index=2*i+p;
	 Norme=Norme+pBem[index]+nBem[index];
	 Normo=Normo+(2.0*i+p+bet)*pBem[index]+(-2.0*i+p+bet)*nBem[index];
	}
     Norme=1.0/Norme;
     Normo=1.0/Normo;
     for(i=0;i<=CDIM;i++)
      {
	 index=2*i+p;
	 pBom[index]=pBem[index]*Normo;
	 nBom[index]=nBem[index]*Normo;
	 pBem[index]=pBem[index]*Norme;
	 nBem[index]=nBem[index]*Norme;
      }
     }
   if(MF==2)
    {
       /*
	   McLachlan's  Normalization convension   p.82
	   The even=odd normalization
       */
     index=p;  // i=0
     if(p==0)
     {
	Norme=pBem[index]*pBem[index];
     }
   else
     {
	 Norme=pBem[index]*pBem[index]+nBem[index]*nBem[index];
     }
    for(i=1;i<=CDIM;i++)
	{
	 index=2*i+p;
	 Norme=Norme+pBem[index]*pBem[index]+nBem[index]*nBem[index];  // McLachlan normalization
	}
     Norme=sqrt(Norme);
     Norme=1.0/Norme;
     Normo=Norme;   // 	   The even=odd normalization
     cr=0;
     for(i=0;i<=CDIM;i++)
      {
	 index=2*i+p;
	 pBom[index]=pBem[index]*Normo;
	 nBom[index]=nBem[index]*Normo;
	 pBem[index]=pBem[index]*Norme;
	 nBem[index]=nBem[index]*Norme;
	 u=u+pBem[index]*pBem[index]+nBem[index]*nBem[index];
	 //  u must be  one  i.e. u=1
      }
     }


   delete[] nW;
   delete[] pW;

   return acc;

}

/*-----------------------------------------------------------------------*/
double CoefficientFr(complexn c, complexn v ,complexn h,
		     complexn pBem[],complexn nBem[],complexn pBom[],complexn nBom[],int CDIM, int MF)
{
    /*
      You can chose several normalization
      0: No normalization
      1: Morse and Feshbach
      2: McLachlan (p.112)

      For coefficient computation
      inputs:
      c=MCN
      v=fractional order
      h=pamarater
      CDIM= Coefficients Array dimension

      outputs:
      pBem[]=Array of Coefficients of positive index for even functions
      nBem[]=Array of Coefficients of negative index for even functions
      pBom[]=Array of Coefficients of positive index for odd  functions
      nBom[]=Array of Coefficients of negative index for odd  functions
      local:
      pW= positive indexed array
      nW= negative indexed array
    */
    complexn *nW,*pW;
    complexn h2,hi2,bet,cr;
    double acc;
    complexn c3,cm1,nC1,pC1;
    int Ml,p,i,vi,r,m;
    int index;
    complexn Norme,Normo;
    complexn Vm,Vm1,Ti,u;

    vi=int(real(v));
    bet=v-complexn(vi,0);
    p=vi%2;


    if(real(c)>=0) cr=0.5*(sqrt(c)-bet-complexn(p,0)); else cr=0;
    r=int(real(cr));  /* Center the computation around the singularity */
    if(r<0) r=0;
    Ml=5+r+CDIM;  // Large m

   h2=0.25*h*h;
   hi2=1.0/h2;
    nW=new complexn [2*CDIM+6];
    pW=new complexn [2*CDIM+6];
    if(nW==NULL)
      {
       cout<<"\n CoefficientFr: out of memory..."; exit(0);
      }
    if(pW==NULL)
      {
       cout<<"\n CoefficientFr: out of memory..."; exit(0);
      }

    for(i=0;i<2*CDIM+2;i++)
     {
       nW[i]=0;
       pW[i]=0;
     }

     /* computing for negative m */
  Vm1=Vm=0;
  for(m=Ml;m>1;m--)  // backward  down to m=2
    {
      u=complexn(-2.0*m+p,0)+bet;
      Ti=1.0/(u*u);
      Vm1=-Ti*h2/(complexn(1,0)-c*Ti+Ti*h2*Vm);
      Vm=Vm1;
      if(m<=CDIM+1) nW[2*m-p-2]=Vm1;
    }
   m=1;              // backward    m=1
   u=complexn(-2.0*m+p,0)+bet;
   Ti=1.0/(u*u);
   Vm1=-Ti*h2/(complexn(1,0)-c*Ti+Ti*h2*Vm);
   nC1=Vm1;

     /* computing for positive m */
   u=bet+complexn(p,0);
   Vm=(c-u*u)/h2-Vm1;  // forward  m=0;
   Vm1=Vm;
   pW[p]=Vm;
  for(m=1;m<=r;m++)  // forward up to m=r
   {
    u=complexn(2.0*m+p,0)+bet;
    Ti=u*u;
    Vm=c*hi2-Ti*hi2-1.0/Vm1;  // Vm1=V_{m-1}  and Vm=V_{m}
    Vm1=Vm;
    if(m<=CDIM) pW[2*m+p]=Vm;
   }

  Vm=0;
  for(m=Ml;m>r;m--)  // backward  down to m=r+1
    {
      u=complexn(2.0*m+p,0)+bet;
      Ti=1.0/(u*u);
      Vm1=-Ti*h2/(complexn(1,0)-c*Ti+Ti*h2*Vm);
      Vm=Vm1;
      if(m<=CDIM) pW[2*m+p-2]=Vm;
    }

      pBem[p]=1;
      for(m=0;m<=CDIM;m++)
       {pBem[2*m+2+p]=pW[2*m+p]*pBem[2*m+p];}

      nBem[0]=0;
      nBem[2-p]=pBem[p]*nC1;
      for(m=1;m<=CDIM;m++)
       {nBem[2*m+2-p]=nW[2*m-p]*nBem[2*m-p];}

	 /* accuracy check see McLachlan p.113 eq(19), at m=0 */
	 /*     c1=pBm[p];c3=pBem[2+p]; cm1=nBem[2-p];  */
     u=bet+complexn(p,0);
     acc=abs((pBem[2+p]+nBem[2-p])*h2-(c-u*u)*pBem[p]);
     if(acc>1e-14) {cout<<"\n Warning: CoefficientFr() is not accurate acc="<<acc;}

   /*----------------------------------------------------------------*/

       /* Normalization convensions */
    if(MF==1)
    {      /*  Morse and Feshbach Normalization convension
	       in MF the even is different from the odd normalization
	   */
     index=p;  // i=0
     if(p==0)
     {
	Norme=pBem[index];
        u=bet+complexn(p,0);
	Normo=u*(pBem[index]);
     }
   else
    {
	 Norme=pBem[index]+nBem[index];  // add C[-1] and C[+1]
         u=bet+complexn(p,0);
	 Normo=u*pBem[index]+u*nBem[index];
    }
     for(i=1;i<=CDIM;i++)
	{
	 index=2*i+p;
	 Norme=Norme+pBem[index]+nBem[index];
	 Normo=Normo+(complexn(2.0*i+p,0)+bet)*pBem[index]+(complexn(-2.0*i+p,0)+bet)*nBem[index];
	}
     Norme=1.0/Norme;
     Normo=1.0/Normo;
     for(i=0;i<=CDIM;i++)
      {
	 index=2*i+p;
	 pBom[index]=pBem[index]*Normo;
	 nBom[index]=nBem[index]*Normo;
	 pBem[index]=pBem[index]*Norme;
	 nBem[index]=nBem[index]*Norme;
      }
     }
   if(MF==2)
    {
       /*
	   McLachlan's  Normalization convension   p.82
	   The even=odd normalization
       */
     index=p;  // i=0
     if(p==0)
     {
	Norme=pBem[index]*pBem[index];
     }
   else
     {
	 Norme=pBem[index]*pBem[index]+nBem[index]*nBem[index];
     }
    for(i=1;i<=CDIM;i++)
	{
	 index=2*i+p;
	 Norme=Norme+pBem[index]*pBem[index]+nBem[index]*nBem[index];  // McLachlan normalization
	}
     Norme=sqrt(Norme);
     Norme=1.0/Norme;
     Normo=Norme;   // 	   The even=odd normalization
     cr=0;
     for(i=0;i<=CDIM;i++)
      {
	 index=2*i+p;
	 pBom[index]=pBem[index]*Normo;
	 nBom[index]=nBem[index]*Normo;
	 pBem[index]=pBem[index]*Norme;
	 nBem[index]=nBem[index]*Norme;
	 u=u+pBem[index]*pBem[index]+nBem[index]*nBem[index];
	 //  u must be  one  i.e. u=1
      }
     }


   delete[] nW;
   delete[] pW;

   return acc;

}

#undef complexn

