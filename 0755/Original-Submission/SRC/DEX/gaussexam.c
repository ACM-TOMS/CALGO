/*
   --------------------------------------------------------------
   File gaussexam.c of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
   gaussexam.c contains the example with gaussian elimination
   described in the manual
*/

#include <stream.h>
#include <stdio.h>
#ifdef __GNUG__
#include <std.h>
#else
#include <stdlib.h>
#endif

#include "adouble.h"    // These includes provide the compiler with
#include "adutils.h"    // definitions and utilities for `adoubles'.

void gausselim(int n, adoublem& A, adoublev& bv)
{
  along i;              // active integer declaration
  adoublev temp(n);     // active vector declaration 
  adouble r,rj,temps;
  int j,k;
  for (k=0; k < n; k++) // elimination loop 
  {
    i = k;
    r = fabs(A[k][k]);  // initial pivot size 
    for (j=k+1; j<n; j++)
    {
      rj = fabs(A[j][k]); 
      // look for a larger element in the same column 
      condassign(i,rj-r,j); // conditional assignment 
      condassign(r,rj-r,rj);
    } // endfor
    temp = A[i]; // switch rows using active subscripting
    A[i] = A[k]; // necessary even if i happens to equal 
    A[k] = temp; // k during taping
    temps = bv[i];
    bv[i]=bv[k];
    bv[k]=temps;
    if (!value(A[k][k])) // passive subscripting
      exit(1); // Matrix singular! 
    temps= A[k][k];
    A[k] /= temps;
    bv[k] /= temps;
    for (j=k+1; j<n; j++)
    {
      temps= A[j][k];
      A[j] -= temps*A[k]; // vector operations
      bv[j] -= temps*bv[k];
    } // endfor
  } // endfor elimination loop
  temp=0.0;
  for(k=n-1; k >= 0; k--)  // backsubstitution
     temp[k] = (bv[k]-(A[k]*temp))/A[k][k];
  bv=temp;
  return;
} // end gausselim

void main() 
{
  int i,j;
  short tag = 1;
  int dum=1;
  const int size=5;
  const int indep=size*size+size;
  const int depen=size;
  double* arguments=new double[indep];
  double* taylors=new double[depen];
  double yp[size],xp[size*size+size];  // passive variable
  double **A_1, **A_2, *a_1, *a_2, *b_1, *b_2;
  A_1=(double**)malloc(size*sizeof(double*));
  A_2=(double**)malloc(size*sizeof(double*));
  a_1=(double*)malloc(size*sizeof(double));
  a_2=(double*)malloc(size*sizeof(double));
  b_1=(double*)malloc(size*sizeof(double));
  b_2=(double*)malloc(size*sizeof(double));
  for(i=0;i<size;i++){
    A_1[i]=(double*)malloc(size*sizeof(double));
    A_2[i]=(double*)malloc(size*sizeof(double));
    a_1[i]=i*0.25;
    a_2[i]=i*0.5;
    b_1[i]=i*0.33;
    b_2[i]=i*0.68;
  } // endfor 
  for(i=0;i<size;i++) {
    for(j=0;j<size;j++) {
      A_1[i][j]=a_1[i]*b_1[j];
      A_2[i][j]=a_2[i]*b_2[j];
    } // endfor 
    A_1[i][i]+=i+1;
    A_2[i][i]+=i+1.5;
  } // endfor
  adoublem A(size,size);
  adoublev bv(size);  // active variables
  int N=size*size;
  trace_on(tag,dum);   // Begin taping all calculations with 'adoubles'
  for(i=0;i<size;i++)
  {
    for(j=0;j<size;j++)
    {
      A[i][j]<<=A_1[i][j]; /* indep. vars */
      arguments[i*size+j]=A_1[i][j]; /* args for forward */
    } /* endfor */
  } /* endfor */
  for(i=0;i<size;i++)
  {
    bv[i]<<=-i-1; /* indep. vars */
    arguments[N+i]=-i-1; /* args for forward */
  }
  gausselim(size,A,bv);
  bv >>= yp;
  trace_off(); 
  forward(tag,depen,indep,0,1,arguments,taylors);
  cout << "Compare the calculated solution components of the\nforward sweep and the direct evaluation: forward - direct = 0 ?\n";
  for(i=0;i<depen;i++)
  cout << taylors[i] << " - " << yp[i] << " = " << taylors[i] - yp[i] << "\n";
  // use the same tape for an other system matrix:
  for(i=0;i<size;i++)
  {
    for(j=0;j<size;j++)
    {
      arguments[j*size+i]=A_2[i][j]; /* new args for forward */
      A[i][j]=A_2[i][j]; /* new args for gausselim */
    } /* endfor */
    bv[i]=-i-1; /* old indep. vars */
  } /* endfor */
  gausselim(size,A,bv); //calculation without taping
  forward(tag,depen,indep,0,1,arguments,taylors);
  cout << "\nThe same comparison for a different system matrix: \n";
  for(i=0;i<depen;i++)
  cout << taylors[i] << " - " << value(bv[i]) << " = " << taylors[i] - value(bv[i]) << "\n";
} // endmain

