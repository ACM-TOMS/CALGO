/*
   --------------------------------------------------------------
   File scalexam.c of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
   scalexam.c contains the scalar example described in the
   manual
*/

#include "adouble.h"
#include "adutils.h"
#include <stream.h>

adouble power(adouble x, int n) {
adouble z=1;
if (n>0) {                         // Recursion and branches
  int nh =n/2;                     // that do not depend on
  z = power(x,nh);                 // adoubles are fine !!!!
  z *= z;
  if (2*nh != n) 
    z *= x;
  return z; 
} // end if
else {
  if (n==0)                        // The local adouble z dies
    return z;                      // as it goes out of scope.
  else 
    return 1/power(x,-n); 
} // end else    
} // end power

void main() {
int i,tag=1;
cout<<"monomial degree=? \n";      // Input the desired degree.
int n; cin >> n;
/*Allocations and Initializations*/
double* Y[1];
*Y = new double[n+2];
double* X[1];                      // Allocate passive variables with
*X = new double[n+4];              // extra dimension for derivatives
X[0][0] = 0.5;                     // function value = 0. coefficient
X[0][1] = 1.0;                     // first derivative = 1. coefficient
for(i=0; i < n+2; i++)
  X[0][i+2]=0;                     // further coefficients.
double* Z[1];                      // used for checking consistency
*Z = new double[n+2];              // between forward and reverse
adouble y,x;                       // Declare active variables
/*Beginning of Active Section*/
trace_on(1);                       // tag = 1 and keep = 0
x <<= X[0][0];                     // Only one independent var
y = power(x,n);                    // Actual function call
y >>= Y[0][0];                     // Only one dependent adouble
trace_off();                       // No global adouble has died
/*End of Active Section */
double u[1];                       // weighting vector
u[0]=1;                            // for reverse call
for(i=0; i < n+2; i++) {           // Note that keep = i+1 in call
  forward(tag,1,1,i,i+1,X,Y);      // Evaluate the i-the derivative
  if (i==0)
    cout << Y[0][i] << " - " << value(y) << " = " << Y[0][i]-value(y)
    << " (should be 0)\n";
  else
    cout << Y[0][i] << " - " << Z[0][i] << " = " << Y[0][i]-Z[0][i] 
    << " (should be 0)\n";
  reverse(tag,1,1,i,u,Z);         // Evaluate the (i+1)-st deriv.
  Z[0][i+1]=Z[0][i]/(i+1);        // Scale derivative to Taylorcoeff.
} // end for 
} // end main

