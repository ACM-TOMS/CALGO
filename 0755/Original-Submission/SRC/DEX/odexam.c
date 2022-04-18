/*
   --------------------------------------------------------------
   File odeexam.c of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
   odeexam.c contains the ODE example described in the manual
*/

#include "adouble.h"
#include "adutils.h"
#include <stream.h>

void tracerhs(short int tag, double* py, double* pyprime) {
adoublev y(3);             //This time we left the parameters
adoublev yprime(3);        // passive and use the vector types.
trace_on(tag);
y <<= py;                  //Initialize and mark independents
yprime[0] = -sin(y[2]) + 1e8*y[2]*(1-1/y[0]);
yprime[1] = -10*y[0] + 3e7*y[2]*(1-y[1]);
yprime[2] = -yprime[0] - yprime[1];
yprime >>= pyprime;        //Mark and pass dependents
trace_off(tag);
} // end tracerhs

void main() {
int i,j,deg;  
int n=3;
double py[3];
double pyp[3];
cout << "degree of Taylor series =?\n";
cin >> deg;
double **X;
X=(double**)malloc(n*sizeof(double*));
for(i=0;i<n;i++)
  X[i]=(double*)malloc((deg+1)*sizeof(double));
double*** Z=new double**[n];
double*** B=new double**[n];
short** nz = new short*[n];
for(i=0;i<n;i++) {
  Z[i]=new double*[n];
  B[i]=new double*[n];
  for(j=0;j<n;j++){
    Z[i][j]=new double[deg];
    B[i][j]=new double[deg];
  } // end for 
} // end for 
for(i=0;i<n;i++) {
  py[i] = (i == 0) ? 1.0 : 0.0;       // Initialize the base point
  X[i][0] = py[i];                    // and the Taylor coefficient;
  nz[i] = new short[n];               // set up sparsity array
} // end for                        
tracerhs(1,py,pyp);                   // trace RHS with tag = 1
forode(1,n,deg,X);                    // compute deg coefficients
reverse(1,n,n,deg-1,Z,nz);            // U defaults to the identity
accode(n,deg-1,Z,B,nz);
cout << "nonzero pattern:\n";
for(i=0;i<n;i++) {
  for(j=0;j<n;j++)
    cout << nz[i][j]<<"\t";
  cout <<"\n"; 
} // end for 
} // end main
