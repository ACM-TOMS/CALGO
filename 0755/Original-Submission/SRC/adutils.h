/*
  ------------------------------------------------------------- 
  File adutils.h of ADOL-C version 1.6 as of January 1,   1995 
  Included in ---> drivers.c

  -------------------------------------------------------------
  Contains the definitions of the standard utility functions for
  ADOL-C.  All of these routines are defined in drivers.c.
  (See adutilsc.h for utility routines that are callable from FORTRAN,
   C, and C++.)
*/

extern "C" {
#include "adutilsc.h"
}
/* Memory management routines. */

extern double ** myalloc(int,int);
extern double *** myalloc(int,int,int);


/* forward(tag,m,n,d,keep,X[n][d+1],Y[m][d+1]) */

void forward(short,int,int,int,int,double**,double**);

/* Y can be one dimensional if m=1 */

void forward(short,int,int,int,int,double**,double*);

/* X and Y can be one dimensional if d = 0 */

void forward(short,int,int,int,int,double*,double*); 

/* forward(tag,m,n,d,p,x[n],X[n][p][d],y[m],Y[m][p][d])        */

void forward(short, int,int,int,int,double*,double***,double*,double***);
 
/* forward(tag,m,n,p,x[n],X[n][p],y[m],Y[m][p])        */

void forward(short, int,int,int,double*,double**,double*,double**);

/* reverse(tag,m,n,d,u[m],Z[n][d+1]) */

void reverse(short,int,int,int,double*,double**);

/* u can be a scalar if m=1; */

void reverse(short,int,int,int,double,double**);

/* Z can be vector if d = 0; Done by specialized code */

void reverse(short,int,int,int,double*,double*);

/* u and Z can be scalars if m=1 and d=0; */

void reverse(short,int,int,int,double,double*);

/* reverse(tag,m,n,d,p,U[m][p],Z[p][n][d+1],nz[p][n]) */

void reverse(short,int,int,int,int,double**,double***,short** =0);

/* U can be a vector if m=1 */
void reverse(short,int,int,int,int,double*,double***,short** = 0);

/* If d=0 then Z may be matrix; Done by specialized code */
void reverse(short,int,int,int,int,double**,double**);

/* If m=1 and d=0 then U can be vector and Z a matrix but no nz. */

void reverse(short,int,int,int,int,double*,double**);

/* If p and U are omitted they default to m and I so that as above */

void reverse(short,int,int,int,double***,short** =0);

/* forode(tag,n,tau,dold,dnew,X[n][d+1]) */

void forode(short,int,double,int,int,double**);

/* the scaling defaults to 1 */

void forode(short,int,/*1.0*/ int,int,double**);

/* previous order defaults to 0 */

void forode(short,int,double,/*0*/ int,double**);

/* both default */

void forode(short,int,/*1.0 , 0*/int,double**);

/* accode(n,tau,d,Z,B,nz) */

void accode(int,double,int,double***,double***,short** = 0 );

/* scaling defaults to 1 */

void accode(int,/*1.0*/ int,double***,double***,short** = 0 );






