/*
   --------------------------------------------------------------
   File vectexam.c of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
   vectexam.c contains the vector example described in the
   manual
*/

#include "adouble.h"
#include "adutils.h"
#include <stream.h>
#include <math.h>

void main() {
int n,i,j,counts[12];
cout << "number of independent variables = ?  \n";
cin >> n;
double* xp = new double[n];          
adouble* x = new adouble[n];         // or: adoublev x(n);
for(i=0;i<n;i++)
  xp[i] = (i+1.0)/(2.0+i);           // some initialization
trace_on(1);                         // tag =1, keep=0 by default
adouble y = 1;
for(i=0;i<n;i++){
  x[i] <<= xp[i];                    // or  x<<= xp outside the loop
  y *= x[i]; 
} // end for
double yp=0.0;
y >>= yp;
delete[] x;                          // Not needed if x adoublev
trace_off();
tapestats(1,counts);                 // Reading of tape statistics
cout<<"maxlive "<<counts[2]<<"\n";
// ..... print other tape stats
double* g = new double[n];           // or: doublev g(n);
gradient(1,n,xp,g);                  // gradient evaluation
double** H=(double**)malloc(n*sizeof(double*));
for(i=0;i<n;i++)
  H[i]=(double*)malloc((i+1)*sizeof(double)); 
hessian(1,n,xp,H);                   // H equals (n-1)g since g is
double errg =0;                      // homogeneous of degree n-1.
double errh =0;
for(i=0;i<n;i++) 
  errg += fabs(g[i]-yp/xp[i]);       // vanishes analytically.
for(i=0;i<n;i++) {
  for(j=0;j<n;j++) {
    if (i>j)                         // lower half of hessian
      errh += fabs(H[i][j]-g[i]/xp[j]); 
  } // end for 
} // end for
cout << yp-1/(1.0+n) << " error in function \n";
cout << errg <<" error in gradient \n";
cout << errh <<" consistency check \n";
} // end main


