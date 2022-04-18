/*
  ------------------------------------------------------------- 
  File adutilsc.h of ADOL-C version 1.6 as of January 1,   1995 
  Included in ---> driversc.c
                   adutils.h


  -------------------------------------------------------------


 
   This header file list the user utilities that are callable from C
   Overloaded C++ routines are listed in adutils.h.  
   Functions prototyped here are defined in the files
            ---> driversc.c
                 hos_forward.c  (hos_forward)
		 fos_reverse.c  (fos_reverse)
		 fov_reverse.c  (fov_reverse)
		 hov_reverse.c  (hov_reverse)
		 hos_reverse.c  (hos_reverse)
		 tapestats      (taputil3.c)

*/


/* tapestats(tag,counts)  from (taputil3.c) */
void tapestats(short,int*); 

typedef long   fint;
typedef double fdouble;

double myclock();
double** myalloc2(int, int);
double*** myalloc3(int, int, int);

/* hos_forward(tag,m,n,d,keep,X[n][d+1],Y[m][d+1]) from (hos_forward.c) */

void hos_forward(short,int,int,int,int,double**,double**);

/* now pack the arrays into vectors for Fortran calling */

fint hos_forward_(fint*,fint*,fint*,fint*,fint*,fdouble*,fdouble*);

/* hov_forward(tag,m,n,d,p,x[n],X[n][p][d],y[m],Y[m][p][d])        */

void hov_forward(short, int,int,int,int,double*,double***,double*,double***);

/* now pack the arrays into vectors for Fortran calling      */

fint hov_forward_(fint*,fint*,fint*,fint*,fint*,fdouble*,fdouble*,fdouble*,fdouble*);

/* fov_forward(tag,m,n,p,x[n],X[n][p],y[m],Y[m][p])   */
 
void fov_forward(short, int,int,int,double*,double**,double*,double**);
 
/* now pack the arrays into vectors for Fortran calling      */

fint fov_forward_(fint*,fint*,fint*,fint*,fdouble*,fdouble*,fdouble*,fdouble*);

/*  hos_reverse(tag,m,n,d,u[m],Z[n][d+1]) from (hos_reverse.c) */

void hos_reverse(short,int,int,int,double*,double**);

/* now pack the arrays into vectors for Fortran calling      */

fint hos_reverse_(fint*,fint*,fint*,fint*,fdouble*,fdouble*);

/* fos_reverse(tag,m,n,u[m],z[n]);  from (fos_reverse.c)     */
void fos_reverse(short,int,int,double*,double*);

/* now pack the arrays into vectors for Fortran calling      */
fint fos_reverse_(fint*,fint*,fint*,fdouble*,fdouble*);

/* hov_reverse(tag,m,n,d,p,U[m][p],Z[p][n][d+1],nz[p][n]) from(hov_reverse.c)*/

void hov_reverse(short,int,int,int,int,double**,double***,short**);

/* now pack the arrays into vectors for Fortran calling      */
fint hov_reverse_(fint*,fint*,fint*,fint*,fint*,fdouble*,fdouble*);

/* fov_reverse(tag,m,n,d,p,U[m][p],Z[p][n]); from (fov_reverse.c)  */

void fov_reverse(short,int,int,int,double**,double**);

/* now pack the arrays into vectors for Fortran calling      */

fint fov_reverse_(fint*,fint*,fint*,fint*,fdouble*,fdouble*);

/* function(tag,m,n,x[n],y[m])                               */

void function(short,int,int,double*,double*);
fint function_(fint*,fint*,fint*,fdouble*,fdouble*);

/* gradient(tag,n,x[n],g[n])                              */

void gradient(short,int,double*,double*);
fint gradient_(fint*,fint*,fdouble*,fdouble*);

/* lagra_hess_vec(tag,m,n,x[n],u[m],v[n],w[n]);              */

void lagra_hess_vec(short,int,int,double*,double*,double*,double*);
fint lagra_hess_vec_(fint*,fint*,fint*,fdouble*,fdouble*,fdouble*,fdouble*);

/* hess_vec(tag,n,x[n],v[n],w[n]);                           */

void hess_vec(short,int,double*,double*,double*);
fint hess_vec_(fint*,fint*,fdouble*,fdouble*,fdouble*);

/* hessian(tag,n,x[n], lower triangle of H[n][n])            */

void hessian(short,int,double*,double**);

/* now pack the arrays into vectors for Fortran calling      */

fint hessian_(fint*,fint*,fdouble*,fdouble*);

/* jacobian(tag,m,n,x[n],J[m][n])                            */

void jacobian(short,int,int,double*,double**);

/* now pack the arrays into vectors for Fortran calling      */

fint jacobian_(fint*,fint*,fint*,fdouble*,fdouble*);

/* jac_vec(tag,m,n,x[n],v[n],u[m]);                          */

void jac_vec(short,int,int,double*,double*,double*);
fint jac_vec_(fint*,fint*,fint*,fdouble*,fdouble*,fdouble*);

/* vec_jac(tag,m,n,repeat,x[n],u[m],v[n]);                   */

void vec_jac(short,int,int,int,double*,double*,double*);
fint vec_jac_(fint*,fint*,fint*,fint*,fdouble*,fdouble*,fdouble*);

/* forodec(tag,n,tau,dold,dnew,X[n][d+1])                    */

void forodec(short,int,double,int,int,double**);
fint forodec_(fint*,fint*,fdouble*,fint*,fint*,fdouble*);

/* accodec(n,tau,d,Z[n][n][d+1],B[n][n][d+1],nz[n][n])       */

void accodec(int,double,int,double***,double***,short**);

/* now pack the arrays into vectors for Fortran calling      */

fint accodec_(fint*,fdouble*,fint*,fdouble*,fdouble*);
 










