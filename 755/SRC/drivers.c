/*
  ------------------------------------------------------------------------
  file drivers.c of ADOL-C version 1.6 as of January 1,   1995  
  ------------------------------------------------------------------------
  This file contains definitions for the functions prototyped in
  adutils.h.  Each function is an ADOL-C C++ utility.
*/


#include "adutils.h"   /* Prototypes */
#include "dvlparms.h"  /* Developers Parameters */


#define fabs(x) ((x) > 0 ? (x) : -(x))
#define ceil(x) ((int)((x)+1) - (int)((x) == (int)(x)))


double** myalloc(int m,int n)
{
  
  /* This function allocates row matrices contiguously  */
  
  return myalloc2(m,n);
  
  /* To deallocate an array set up by   A = myalloc(m,n)  */
  
  /*  use  free((char*)*A); free((char*)A); in that order */
  
}

double*** myalloc(int m,int n,int p)
{
  
  /* This function allocates 3-tensors contiguously */
  
  return  myalloc3(m,n,p);
  
  /* To deallocate an array set up by  A = myalloc(m,n,p)     */

  /* use free((char*)**A); free((char*)*A); free((char*)A) ;  */
  
}

void forward(short tag,
             int depen,
             int indep,
             int degre,
             int keep,
             double** X,
             double** Y)
{
  hos_forward(tag,depen,indep,degre,keep,X,Y);
}

void forward(short tag,
	     int depen,
	     int indep,
	     int degre,
	     int keep,
	     double** X,
	     double* Y)
{
  if(depen==1) 
    hos_forward(tag,depen,indep,degre,keep,X,&Y);
  else 
    {
      printf("ADOL-C error: wrong Y dimension in forward \n");
      exit(-1);
    }
}

void forward(short tag,
	     int depen,
	     int indep,
	     int degre,
	     int keep,
	     double* X,
	     double* Y)
{
  static double **Xl, **Yl;
  static int indax, depax;
  if(degre != 0)
    {
      printf("ADOL-C error:  wrong X and Y dimensions in forward \n");
      exit(-1);
    }
  else
    {
      if (indep compsize indax || depen compsize depax)
	{
	  if(indax)
	    {
	      free((char*)*Xl);
	      free((char*)Xl);
	      free((char*)*Yl);
	      free((char*)Yl);
	    }
	  Xl = myalloc(indep,1);
	  Yl = myalloc(depen,1);
	  indax = indep;
	  depax = depen;
	}
    }
  for(int i=0; i<indep; i++)
    *Xl[i] = X[i];
  hos_forward(tag,depen,indep,degre,keep,Xl,Yl);
  for(i=0; i<depen; i++)
    Y[i] = *Yl[i];
}

void forward(short tag,
             int depen,
             int indep,
             int degre,
             int taylnum,
             double* base,
             double*** X,
             double* value,
             double*** Y)
{
  hov_forward(tag,depen,indep,degre,taylnum,base,X,value,Y);
}

void forward(short tag,
             int depen,
             int indep,
             int taylnum,
             double* base,
             double** X,
             double* value,
             double** Y)
{
  fov_forward(tag,depen,indep,taylnum,base,X,value,Y);
}


void reverse(short tag, 
             int depen,
	     int indep,
	     int degre,
	     double* u,
	     double** Z)
{
  hos_reverse(tag,depen,indep,degre,u,Z);
}

void reverse(short tag, 
             int depen,
	     int indep,
	     int degre,
	     double u,
	     double** Z)
{
  if(depen != 1) 
    {
      printf("ADOL-C error:  wrong u dimension in scalar-reverse \n");
      exit(-1);
    }
  else
    {
      hos_reverse(tag,depen,indep,degre,&u,Z);
    }
}

void reverse(short tag, 
             int depen,
	     int indep,
	     int degre,
	     double* u,
	     double* Z)
{
  if(degre != 1) 
    {
      printf("ADOL-C error:  wrong Z dimension in scalar-reverse \n");
      exit(-1);
    }
  fos_reverse(tag,depen,indep,u,Z);
}

void reverse(short tag, 
             int depen,
	     int indep,
	     int degre,
	     double u,
	     double* Z)
{
  if(depen != 1 || degre !=0 ) 
    {
      printf("ADOL-C error:  wrong u or Z dimension in scalar-reverse \n");
      exit(-1);
    }
  else
    fos_reverse(tag,depen,indep,&u,Z);
}

void reverse(short tag, 
             int depen,
	     int indep,
	     int degre,
	     int nrows,
	     double** U,
	     double*** Z,
	     short** nz)
{
  hov_reverse(tag,depen,indep,degre,nrows,U,Z,nz);
}

void reverse(short tag, 
             int depen,
	     int indep,
	     int degre,
	     int nrows,
	     double* U,
	     double*** Z,
	     short** nz)
{
  if(depen != 1 ) 
    {
      printf("ADOL-C error:  wrong U dimension in vector-reverse \n");
      exit(-1);
    }
  else
    hov_reverse(tag,depen,indep,degre,nrows,&U,Z,nz);
}

void reverse(short tag, 
             int depen,
	     int indep,
	     int degre,
	     int nrows,
	     double** U,
	     double** Z)
{
  if(degre != 0) 
    {
      printf("ADOL-C error:  wrong degree in vector-reverse \n");
      exit(-1);
    }
  else
    fov_reverse(tag,depen,indep,nrows,U,Z);
}

void reverse(short tag, 
             int depen,
	     int indep,
	     int degre,
	     int nrows,
	     double* U,
	     double** Z)
{
  if(depen != 1) 
    {
      printf("ADOL-C error:  wrong U dimension in vector-reverse \n");
      exit(-1);
    }
  else
    fov_reverse(tag,depen,indep,nrows,&U,Z);
}

void reverse(short tag, 
             int depen,
	     int indep,
	     int degre,
	     double*** Z,
	     short** nz)
{
  static int depax;
  static double** I;
  if(depen compsize depax)
    {
      if(depax) 
	{
	  free((char*)*I); free((char*)I);
	}
      I = myalloc(depen,depen);
      for (int i=0;i<depen;i++)
	{
	  for (int j=0;j<depen;j++)
	    I[i][j] = 0;
	  I[i][i] =1;
	}
      depax = depen;
    }
  hov_reverse(tag,depen,indep,degre,depen,I,Z,nz);
}

void forode(short tag,             // tape identifier
            int n,                 // space dimension
            double tau,            // scaling defaults to 1.0
            int dol,               // previous degree defaults to zero
            int deg,               // New degree of consistency
            double** y)            // Taylor series
{
  forodec(tag,n,tau,dol,deg,y);
}


void forode(short tag, int n, double tau, int deg, double** y)
{
  /***   Default for previous degree is zero, do things from scratch ***/
  int zero = 0;
  forodec(tag, n, tau, zero, deg, y);
}

void forode(short tag, int n, int dol, int deg, double** y)
{
  /***   Default for scaling is 1.0       *****/
  double tau = 1.0;
  forodec(tag, n, tau, dol, deg, y);
}

void forode(short tag, int n, int deg, double** y)
{
  /***    Combination of both previous defaults   ******/
  double tau = 1.0;
  forode(tag, n, tau, deg, y);
}

void accode(int n,             // space dimension
            double tau,        // scaling defaults to 1.0
            int deg,           // highest degree
	    double*** A,       // input tensor of "partial" Jacobians
            double*** B,       // output tensor of "total" Jacobians
	    short** nonzero )  // optional sparsity characterization 
{
  accodec(n,tau,deg,A,B,nonzero);
}

void accode(int n,             // space dimension
            int deg,           // highest degree
	    double*** A,       // input tensor of "partial" Jacobians
            double*** B,       // output tensor of "total" Jacobians
	    short** nonzero )  // optional sparsity characterization 
{
  double tau = 1.0;
  accodec(n,tau,deg,A,B,nonzero);
}

