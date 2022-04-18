/*
  ------------------------------------------------------------------------
  file driversc.c of ADOL-C version 1.6 as of January 1,   1995  
  ------------------------------------------------------------------------
  This file contains definitions for the functions prototyped in
  adutilsc.h.  Each function is an ADOL-C straight C utility.
*/

/*--------------------------------------------------*/
/* Included timing data for myclock                 */
/*--------------------------------------------------*/


#ifdef __STDC__
#include <time.h>

double myclock(void)
{ return (double)clock() / CLOCKS_PER_SEC; }

#else
#include <sys/types.h>
#include <sys/times.h>

#ifndef HZ	/* clock ticks per second */
#define HZ 60.
#endif

double myclock() {
  struct tms t;
  times(&t);
  return t.tms_utime / HZ;
} /* end myclock */
#endif


#ifdef __cplusplus
extern "C" {
#endif

#include "dvlparms.h" /* Developers Parameters */
#include "adutilsc.h" /* Function Prototypes   */



double* myalloc1(int m)
{
  double* A = (double*)malloc(m*sizeof(double));
  return A;
}

double** myalloc2(int m,int n)
{
  double* Adum = (double*)malloc(m*n*sizeof(double));
  double**  A = (double**)malloc(m*sizeof(double*));
  int i;
  for(i=0;i<m;i++)
    {
      A[i] = Adum;
      Adum += n;
    }
  return A;

  /* To deallocate an array set up by   A = myalloc2(m,n)   */
  /*   use  free((char*)*A); free((char*)A); in that order  */

 }

double*** myalloc3(int m,int n,int p)
{
  /* This function allocates 3-tensors contiguously */ 
  double* Adum = (double*) malloc(m*n*p*sizeof(double));
  double**   Apt = (double**)malloc(m*n*sizeof(double*));
  double***   A = (double***)malloc(m*sizeof(double**));
  int i,j;
  for(i=0;i<m;i++)
    {
      A[i] = Apt;
      for(j=0;j<n;j++)
	{
	  *Apt++ =  Adum;
	  Adum += p;
	}
    }   
  return A;

 /* To deallocate an array set up by  A = myalloc3(m,n,p)  */

  /* use free((char*)**A); free((char*)*A); free((char*)A); */

 }

static void spread1(int m,
		    fdouble* x,
		    double* X)
{ 
  int j;
  for(j=0;j<m;j++)
    X[j] = *x++;
}

static void pack1(int m,
		  double* X,
		  fdouble* x)
{ 
  int j;
  for(j=0;j<m;j++)
      *x++ = X[j] ;
}

static void spread2(int m,
		    int n,
		    fdouble* x,
		    double** X)
{ 
  int i,j;
  for(j=0;j<n;j++)
    for(i=0;i<m;i++)
      X[i][j] = *x++;
}

static void pack2(int m,
		  int n,
		  double** X,
		  fdouble* x)
{ 
  int i,j;
  for(j=0;j<n;j++)
    for(i=0;i<m;i++)
      *x++ = X[i][j] ;
}

static void spread3(int m,
		    int n,
		    int p,
		    fdouble* x,
		    double*** X)
{ 
  int i,j,k;
  for(k=0;k<p;k++)
    for(j=0;j<n;j++)
      for(i=0;i<m;i++)
	X[i][j][k] = *x++;
}

static void pack3(int m,
		  int n,
		  int p,
		  double*** X,
		  fdouble* x)
{
  int i,j,k;
  for(k=0;k<p;k++)
    for(j=0;j<n;j++)
      for(i=0;i<m;i++)
	*x++ = X[i][j][k];
}

fint hos_forward_(fint* ftag,
		 fint* fm,
		 fint* fn,
		 fint* fd,
		 fint* fk,
		 fdouble* fx,
		 fdouble* fy)
{
  int tag=*ftag, m=*fm, n=*fn, d=*fd, k=*fk;
  double** X = myalloc2(n,d+1);
  double** Y = myalloc2(m,d+1);
  spread2(n,d+1,fx,X);
  hos_forward(tag,m,n,d,k,X,Y);
  pack2(m,d+1,Y,fy);
  free((char*)*X);free((char*)X);
  free((char*)*Y);free((char*)Y);
  return 1;
}

fint hov_forward_(fint* ftag,
                 fint* fm,
                 fint* fn,
                 fint* fd,
                 fint* fp,
                 fdouble* fbase,
                 fdouble* fx,
                 fdouble* fvalue,
                 fdouble* fy)
{
  int tag=*ftag, m=*fm, n=*fn, d=*fd, p=*fp;
  double* base = myalloc1(n);
  double* value = myalloc1(m);
  double*** X = myalloc3(n,p,d);
  double*** Y = myalloc3(m,p,d);
  spread1(n,fbase,base);
  spread3(n,p,d,fx,X);
  hov_forward(tag,m,n,d,p,base,X,value,Y);
  pack3(m,p,d,Y,fy);
  pack1(m,value,fvalue);
  free((char*)**X); free((char*)*X); free((char*)X);
  free((char*)**Y); free((char*)*Y); free((char*)Y);
  free((char*)base); free((char*)value);
  return 11;
}

fint fov_forward_(fint* ftag,
                 fint* fm,
                 fint* fn,
                 fint* fp,
                 fdouble* fbase,
                 fdouble* fx,
                 fdouble* fvalue,
                 fdouble* fy)
{
  int tag=*ftag, m=*fm, n=*fn, p=*fp;
  double* base = myalloc1(n);
  double* value = myalloc1(m);
  double** X = myalloc2(n,p);
  double** Y = myalloc2(m,p);
  spread1(n,fbase,base);
  spread2(n,p,fx,X);
  fov_forward(tag,m,n,p,base,X,value,Y);
  pack2(m,p,Y,fy);
  pack1(m,value,fvalue);
  free((char*)*X); free((char*)X);
  free((char*)*Y); free((char*)Y);
  free((char*)base); free((char*)value);
  return 12;
}

  

fint hos_reverse_(fint* ftag,
		 fint* fm,
		 fint* fn,
		 fint* fd,
		 fdouble* fu,
		 fdouble* fz)
{
  int tag=*ftag, m=*fm, n=*fn, d=*fd;
  double** Z = myalloc2(n,d+1);
  double* u = myalloc1(n);
  spread1(n,fu,u);
  hos_reverse(tag,m,n,d,u,Z);
  pack2(n,d+1,Z,fz);
  free((char*)*Z); free((char*)Z);
  free((char*)u);
  return 2;
}

fint fos_reverse_(fint* ftag,
		 fint* fm,
		 fint* fn,
		 fdouble* fu,
		 fdouble* fz)
{
  int tag=*ftag, m=*fm, n=*fn;
  double* u = myalloc1(n);
  double* Z = myalloc1(n);
  spread1(n,fu,u);
  fos_reverse(tag,m,n,u,Z);
  pack1(n,Z,fz);
  free((char*)Z); free((char*)u);
  return 3;
}

fint hov_reverse_(fint* ftag,
		 fint* fm,
		 fint* fn,
		 fint* fd,
		 fint* fp,
		 fdouble* fu,
		 fdouble* fz)
{
  int tag=*ftag, m=*fm, n=*fn, d=*fd, p=*fp;
  double** U = myalloc2(p,m);
  double*** Z = myalloc3(p,n,d+1);
  short ** nop = 0;
  spread2(p,m,fu,U);
  hov_reverse(tag,m,n,d,p,U,Z,nop);
  pack3(p,n,d+1,Z,fz);
  free((char*)**Z); free((char*)*Z); free((char*)Z);
  free((char*)*U); free((char*)U);
  return 4;
}

fint fov_reverse_(fint* ftag,
		 fint* fm,
		 fint* fn,
		 fint* fp,
		 fdouble* fu,
		 fdouble* fz)
{
  int tag=*ftag, m=*fm, n=*fn, p=*fp;
  double** U = myalloc2(p,m);
  double** Z = myalloc2(p,n);
  spread2(p,m,fu,U);
  fov_reverse(tag,m,n,p,U,Z);
  pack2(p,n,Z,fz);
  free((char*)*Z); free((char*)Z);
  free((char*)*U); free((char*)U);
  return 5;
}

void lagra_hess_vec(short tag,
		    int m,
		    int n,
		    double *argument,
		    double *tangent,
		    double *lagrange,
		    double *result)
{ 
  int i;
  int degree = 1;
  int keep = degree+1;
  static double **X, **Y;
  static int maxn, maxm;
  if(n > maxn || m > maxm)
    {
      if(maxn)
	{
	  free((char*)*X); free((char*)X);
	  free((char*)*Y); free((char*)Y);
	}
      X = myalloc2(n,2);
      Y = myalloc2(m,2);
      maxn = n;
      maxm = m;
    } 
  for(i=0;i<n;i++)
    {
      X[i][0] = argument[i];
      X[i][1] = tangent[i];
    }
  hos_forward(tag,m,n,degree,keep,X,Y);
  hos_reverse(tag,m,n,degree,lagrange,X);
  for(i=0;i<n;i++)
    result[i] = X[i][1];
}

fint lagra_hess_vec_(fint* ftag,
		    fint* fm,
		    fint* fn,
		    fdouble *fargument,
		    fdouble *ftangent,
		    fdouble *flagrange,
		    fdouble *fresult)
{
  int tag=*ftag, m=*fm, n=*fn;
  double *argument = myalloc1(n);
  double *tangent = myalloc1(n);
  double *lagrange = myalloc1(m);
  double *result = myalloc1(n);
  spread1(n,fargument,argument);
  spread1(n,ftangent,tangent);
  spread1(m,flagrange,lagrange);
  lagra_hess_vec(tag,m,n,argument,tangent,lagrange,result);
  pack1(n,result,fresult);
  free((char*)argument); free((char*)tangent); free((char*)lagrange);
  free((char*)result);
  return 16;
}

static double one = 1.0;

void hess_vec(short tag,
	      int n,
	      double *argument,
	      double *tangent,
	      double *result)
{
  lagra_hess_vec(tag,1,n,argument,tangent,&one,result);
}

fint hess_vec_(fint* ftag,
	      fint* fn,
	      fdouble *fargument,
	      fdouble *ftangent,
	      fdouble *fresult)
{ 
  int tag=*ftag, n=*fn;
  double *argument = myalloc1(n);
  double *tangent = myalloc1(n);
  double *result = myalloc1(n);
  spread1(n,fargument,argument);
  spread1(n,ftangent,tangent);
  hess_vec(tag,n,argument,tangent,result);
  pack1(n,result,fresult);
  free((char*)argument); free((char*)tangent); free((char*)result);
  return 15;
}

void hessian(short tag,
	     int n,
	     double* argument,
	     double** hess)
{
  int i,j;
  double *v = (double*)malloc(sizeof(double)*n); 
  double *w = (double*)malloc(sizeof(double)*n); 
  for(i=0;i<n;i++)
    v[i] = 0;
  for(i=0;i<n;i++)
    {
      v[i] = 1;   
      hess_vec(tag,n,argument,v,w);
      for(j=0;j<=i;j++)
	hess[i][j] = w[j];
      v[i] = 0;
    }

  free((char *)v); 
  free((char *) w); 
  /* Note that only the lower triangle of hess is filled in */
}

fint hessian_(fint* ftag,
	     fint* fn,
	     fdouble* fx,
	     fdouble* fh) /* length of h should be n*n but the 
                           upper half of this matrix remains unchanged */
{
  int tag=*ftag, n=*fn;
  int i;
  double** H = myalloc2(n,n);
  double* x = myalloc1(n);
  spread1(n,fx,x);
  hessian(tag,n,x,H);
  pack2(n,n,H,fh);
  free((char*)*H); free((char*)H);
  free((char*)x);
  return 7;
}

void jacobian(short tag,
	      int depen,
	      int indep,
	      double *argument,
	      double **jacobian)
{
  int i,j;
  static int mmax, nmax;
  static double **result, **I, **X; 
  if(depen > mmax || indep > nmax)   
    {
      if(mmax) {free((char*)*I);
		free((char*)I);
		free((char*)*result);
		free((char*)result);
		free((char*)*X);
		free((char*)X);}
      I = myalloc2(depen,depen);
      for (i=0;i<depen;i++)
	{
	  for (j=0;j<depen;j++)
	    I[i][j] =0;
	  I[i][i] =1;
	}
      result = myalloc2(depen,1);
      X  = myalloc2(indep,1);
      mmax = depen;
      nmax = indep;
    }
  for(i=0;i<indep;i++)
    *X[i] = argument[i];
  hos_forward(tag,depen,indep,0,1,X,result);
  fov_reverse(tag,depen,indep,depen,I,jacobian);
}

fint jacobian_(fint* ftag,
	      fint* fdepen,
	      fint* findep,
	      fdouble *fargument,
	      fdouble *fjac)
{
  int tag=*ftag, depen=*fdepen, indep=*findep;
  double** Jac = myalloc2(depen,indep);
  double* argument = myalloc1(indep);
  spread1(indep,fargument,argument);
  jacobian(tag,depen,indep,argument,Jac);
  pack2(depen,indep,Jac,fjac);
  free((char*)*Jac); free((char*)Jac);
  free((char*)argument);
  return 8;
}

void jac_vec(short tag,
	     int m, 
	     int n,
	     double* argument,
	     double* tangent,
	     double* column)
{
  int i;
  static double **X, **Y;
  static int maxn, maxm;
  if(n > maxn || m > maxm)
    {
      if(maxn)
	{
	  free((char*)*X); free((char*)X);
	  free((char*)*Y); free((char*)Y);
	}
      X = myalloc2(n,2);
      Y = myalloc2(m,2);
      maxn = n;
      maxm = m;
    } 
  for(i=0;i<n;i++)
    {
      X[i][0] = argument[i];
      X[i][1] = tangent[i];
    }
  hos_forward(tag,m,n,1,0,X,Y);
  for(i=0;i<m;i++)
    column[i] = Y[i][1];
}

fint jac_vec_(fint* ftag,
	     fint* fm, 
	     fint* fn,
	     fdouble* fargument,
	     fdouble* ftangent,
	     fdouble* fcolumn)
{
  int tag=*ftag, m=*fm, n=*fn;
  double* argument = myalloc1(n);
  double* tangent = myalloc1(n);
  double* column = myalloc1(m);
  spread1(n,ftangent,tangent);
  spread1(n,fargument,argument);
  jac_vec(tag,m,n,argument,tangent,column);
  pack1(m,column,fcolumn);
  free((char*)argument); free((char*)tangent); free((char*)column);
  return 17;
}

void vec_jac(short tag,
	     int m,
	     int n,
	     int repeat,
	     double* argument,
	     double* lagrange,
	     double* row)
{
  int i;
  static double **X, **Y;
  static int maxn, maxm;
  if(n > maxn || m > maxm)
    {
      if(maxn)
	{
	  free((char*)*X); free((char*)X);
	  free((char*)*Y); free((char*)Y);
	}
      X = myalloc2(n,1);
      Y = myalloc2(m,1);
      maxn = n;
      maxm = m;
    } 
  for(i=0;i<n;i++)
    {
    X[i][0] = argument[i];
    }
  if(!repeat) hos_forward(tag,m,n,0,1,X,Y);
  fos_reverse(tag,m,n,lagrange,row); 
}

fint vec_jac_(fint* ftag,
	     fint* fm,
	     fint* fn,
	     fint* frepeat,
	     fdouble* fargument,
	     fdouble* flagrange,
	     fdouble* frow)
{ 
  int tag=*ftag, m=*fm, n=*fn, repeat=*frepeat;
  double* argument = myalloc1(n);
  double* lagrange = myalloc1(m);
  double* row = myalloc1(n);
  spread1(m,flagrange,lagrange);
  spread1(n,fargument,argument);
  vec_jac(tag,m,n,repeat,argument,lagrange, row);
  pack1(n,row,frow);
  free((char*)argument); free((char*)lagrange); free((char*)row);
  return 18;
}

void function(short tag,
              int m,
	      int n,
	      double* argument,
	      double* result)
{
  int i,j;
  static double **X, **Y;
  static int maxn, maxm;
  if(n > maxn || m > maxm)
    {
      if(maxn)
	{
	  free((char*)*X); free((char*)X);
	  free((char*)*Y); free((char*)Y);
	}
      X = myalloc2(n,1);
      Y = myalloc2(m,1);
      maxn = n;
      maxm = m;
    } 
  for(i=0;i<n;i++) 
    *X[i] = argument[i];
  hos_forward(tag,m,n,0,0,X,Y);
  for(j=0;j<m;j++)
    result[j] = *Y[j];
  /* free((char*)*X);free((char*)X);
     free((char*)*Y);free((char*)Y); */
}
  
fint function_(fint* ftag,
              fint* fm,
	      fint* fn,
	      fdouble* fargument,
	      fdouble* fresult)
{
  int tag=*ftag, m=*fm,  n=*fn;
  double* argument = myalloc1(m);
  double* result = myalloc1(n);
  spread1(m,fargument,argument);
  function(tag,m,n,argument,result);
  pack1(n,result,fresult);
  free((char*)argument); free((char*)result); 
  return 19;
}

void gradient(short tag,
	      int n,
	      double* argument,
	      double* result)
{
  int i;
  double y2;
  double* y1=&y2;
  double** y=&y1;
  static double **X;
  static int maxn;
  if(n > maxn)
    {
      if(maxn)
	{
	  free((char*)*X); free((char*)X);
	}
      X = myalloc2(n,1);
      maxn = n;
    } 
  for(i=0;i<n;i++)
    *X[i] = argument[i];
  hos_forward(tag,1,n,0,1,X,y);
  vec_jac(tag,1,n,1,argument,&one,result);
}

fint gradient_(fint* ftag,
	      fint* fn,
	      fdouble* fargument,
	      fdouble* fresult)
{
  int tag=*ftag, n=*fn;
  double* argument=myalloc1(n);
  double* result=myalloc1(n);
  spread1(n,fargument,argument);
  gradient(tag,n,argument,result);
  pack1(n,result,fresult);
  free((char*)result); free((char*)argument);
  return 20;
}

void forodec(short tag,    /* tape identifier */
	     int n,        /* space dimension */
	     double tau,   /* scaling defaults to 1.0 */
	     int dol,      /* previous degree defaults to zero */
	     int deg,      /* New degree of consistency        */
	     double** y)   /* Taylor series */
{
  /*********************************************************************
    This is assumed to be the autonomous case.
    Here we are just going around computing the vectors 
    y[][i] for  dol < i <= deg
    by successive calls to forward that works on the tape identified
    by tag. This tape (array of file) must obviously have been
    generated by a the execution of an active section between
    trace_on and trace_off with n independent and n dependent variables
    y must have been set up as  pointer to an array of n pointers
    to double arrays containing at least deg+1 components.
    The scaling by tau is sometimes necessary to avoid overflow.
    **********************************************************************/
  
  static int i, j, k, nax, dax;
  static double **z, taut;
  if ( n > nax || deg > dax )
    {
      if(nax) {free((char*) *z); free((char*) z);}
      z = myalloc2(n,deg+1);     
      nax = n;
      dax = deg;
    }
  
  /******  Here we get  going    ********/
  for (j=dol;j<deg;j++)
    {
      k = (deg)*(j == deg-1 ) ;         /* keep death values in prepration */
      hos_forward(tag,n,n,j,k,y,z);     /* for  reverse called by jacode   */
      taut = tau/(1+j);                 /* only the last time through.     */
      for (i=0;i<n;i++)
	y[i][j+1] = taut*z[i][j];
    }
}

fint forodec_(fint* ftag,    /* tape identifier */
            fint* fn,       /* space dimension */
            fdouble* ftau,  /* scaling defaults to 1.0 */
            fint* fdol,     /* previous degree defaults to zero */
            fint* fdeg,     /* New degree of consistency        */
	    fdouble* fy)   /* Taylor series                    */
{
  int tag=*ftag, n=*fn, dol=*fdol, deg=*fdeg;
  int i;
  double tau=*ftau;
  double** Y = myalloc2(n,deg+1);
  for(i=0;i<n;i++)
    *Y[i] = fy[i];
  forodec(tag,n,tau,dol,deg,Y);
  pack2(n,deg+1,Y,fy);
  free((char*)*Y); free((char*)Y);
  return 9;
}

void accodec(int n,             /* space dimension */
            double tau,         /* scaling defaults to 1.0 */
            int deg,            /* highest degree          */
	    double*** A,        /* input tensor of "partial" Jacobians */
            double*** B,        /* output tensor of "total" Jacobians  */
	    short** nonzero )   /* optional sparsity characterization  */
{
/*
    The purpose of this subroutine is to compute the total derivatives
               B[i=0...n-1][j=0...n-1][k=0...deg].
    The matrix obtained for fixed k represents the Jacobian of the
    (k+1)-st Taylor coefficient vector with respect to the base point of
    the ODE, i.e., the 0-th coefficient vector. The input array
	       A[i=0...n-1][j=0...n-1][k=0...deg]
    has exactly the same format, except that it-s k-th matrix slice
    represents a partial derivative in that the indirect dependence
    of the k-th coefficient vector on the base point via the (k-1)-st
    and other lower Taylor coeffcients has not been taken into account.
    The B's are compute from the A's by the chainrule with the parameter
    tau thrown in for scaling. The calculation is performed so that
    A may directly be overwritten by B i.e. their pointers arguments may
    coincide to save storage.
	   Sparsity is used so far only to reduce the operations count 
    but not to save space. In general we expect that for each given pair
    (i,j) the entries A[i][j][k=0...] are nonzero either for all k, or for no k,
    or for k=0 only. 
	   On entry the optional short array nonzero may be used to identify
    all entries of the A[.][.][k] that are potentially nonzero, i.e. 

            nonzero[i][j] <= 0  implies   A[i][j][k] = 0 for all k 
	    nonzero[i][j] = 1   implies   A[i][j][k] = 0 for all k > 0 .

    In other words we only allow the sparsity of the matrices A[.][.][k]
    to be increasing in that A[.][.][1] is possibly sparser than A[.][.][0]
    and all subseqent A[.][.][k] with k > 0 have the same sparsity pattern.
    That is the typical situation since A[.][.][k] is the k-th
    Taylor coefficient in the time expansion of the Jacobian of the
    right hand side. The entries of this square matrix tend to be either
    constant or trancendental functions of time.
	   The matrices B_k = B[.][.][k] are obtained from the A_k = A[.][.][k]
    by the recurrence
                        tau   /        k                  \
	         B_k = ----- |  A_k + SUM A_{j-1} B_{k-j}  |
	                k+1   \       j=1                 /

    Assuming that the diagonal entries A[i][i][0] are structurally nonzero
    we find that the matrices B[.][.][k=1..] can only lose sparsity
    as k increase. Therfore, we can redefine the nonpositive values 
    nonzero[i][j] so that on exit

	    k <= -nonzero[i][j]  implies     B[i][j][k] = 0 

    which is trivially satisfied for all positive values of nonzero[i][j].
    Due to the increasing sparsity of the A_i and the decreasing sparsity
    of the B_i the first product in the sum of the RHS above determines the
    sparsity pattern of the resulting B_k. Hence the optimal values of
    the nonzero[i][j] depend only on the sparsity pattern of A_0. More
    specifically, all positive -nonzero[i][j] represent the length of the
    shortest directed path connecting nodes j and i in the incidence graph 
    of A_0. 
*/    

  int i,j,k,m,p,nzip,nzpj,isum;
  double *Aip, *Bpj, scale, sum;
  for (k=0;k<=deg;k++)   /* Lets calculate B_k */
    {
      scale = tau/(1.0+k);
      if(nonzero)
	{
	  for (i=0;i<n;i++)
	    for(j=0;j<n;j++)
	      if(k   < -nonzero[i][j]) 
		B[i][j][k] = 0.0;
	      else   
		{
		  sum = A[i][j][k];
		  isum = (nonzero[i][j] >  0);
		  for (p=0;p<n;p++)
		    {
		      nzip = nonzero[i][p];
		      nzpj = nonzero[p][j];
		      if(nzpj > 0) nzpj = 0;
		      if(nzip > 0 && k > -nzpj )  /*otherwise all terms vanish*/
			{
			  Aip = A[i][p];
			  Bpj = B[p][j]+k-1;
			  sum += *Aip*(*Bpj);
			  isum =1;
			  if(nzip > 1 )   /* the A[i][p][m>0] may be nonzero*/
			    for(m=k-1; m>-nzpj;m--)
			      sum += *(++Aip)*(*(--Bpj));
			}
		    }
		  if(isum)         /* we found something nonzero after all*/
		    B[i][j][k] = sum*scale;
		  else
		    {B[i][j][k]= 0;
		     nonzero[i][j]--;
		   }
		}
	}
      else
	{
	  for (i=0;i<n;i++)
	    for(j=0;j<n;j++)
	      {
		sum = A[i][j][k];
		for (p=0;p<n;p++)
		  {
		    Aip = A[i][p];
		    Bpj = B[p][j]+k-1;
		    for(m=k; m>0 ;m--)
		      sum += *(Aip++)*(*Bpj--);
		    B[i][j][k] = sum*scale;
		  }
	      }
	}
    }
}

fint accodec_(fint* fn,             /* space dimension */
	      fdouble* ftau,        /* scaling defaults to 1.0 */
              fint* fdeg,           /* highest degree          */ 
	      fdouble* fa,          /* input tensor of "partial" Jacobians */
              fdouble* fb)          /* output tensor of "total" Jacobians  */
{ 
  int n=*fn, deg=*fdeg;
  double tau=*ftau;
  double*** A = myalloc3(n,n,deg); 
  double*** B = myalloc3(n,n,deg); 
  spread3(n,n,deg,fa,A);
  accodec(n,tau,deg,A,B,0);
  pack3(n,n,deg,B,fb);
  free((char*)**A); free((char*)*A); free((char*)A);
  free((char*)**B); free((char*)*B); free((char*)B);
  return 10;
}

#ifdef __cplusplus
}
#endif

