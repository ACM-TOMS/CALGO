/*
   --------------------------------------------------------------
   File detexam.c of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
   detexam.c contains the determinant example described in the
   manual
*/

#include "adouble.h"
#include "adutils.h"
#include <stream.h>

adouble** A;                        // A is an n x n matrix
int n;                              // k <= n is the order
adouble det(int k, int m) {         // of the submatrix
  if(m == 0 ) return 1.0 ;          // its column indices
  else {                            // are encoded in m.
    adouble* pt = A[k-1];
    adouble t =0 ;
    int s, p =1;
    if (k%2) s = 1; else s = -1;
    for(int i=0;i<n;i++) {
      int p1 = 2*p;
      if ( m%p1 >= p ) {
        t += *pt*s*det(k-1, m-p);   // Recursive call to det.
        s = -s; }
      ++pt;
      p = p1; }
    return t; }
}

void main() {
  int i, m=1,tag=1,keep=1;
  cout << "order of matrix = ? \n"; // Select matrix size
  cin >> n;
  A = new adouble*[n];              // or adoublem A(n,n);
  trace_on(tag,keep);               //      tag=1=keep
  double detout=0.0 , diag = 1.0;   // here keep the intermediates for
  for (i=0; i<n; i++) {             // the subsequent call to reverse
    m *=2;
    A[i] = new adouble[n];          // not needed for adoublem
    adouble* pt = A[i];
    for (int j=0;j<n; j++)
      A[i][j] <<= j/(1.0+i);        //make all elements of A independent
    diag += value(A[i][i]);         //value(adouble) converts to double
    A[i][i] += 1.0; }
  det(n,m-1) >>= detout;            // Actual function call.
  printf("\n %f - %f = %f  (should be 0)\n",detout,diag,detout-diag);
  trace_off();
  double u[1];
  u[0] = 1.0;
  double* B = new double[n*n];
  reverse(tag,1,n*n,1,u,B);
  cout <<" \n first base? : ";
  for (i=0;i<n;i++) {
    adouble sum = 0;
    for (int j=0;j<n;j++)           // The matrix A times the first n
      sum += A[i][j]*B[j];          // components of the gradient B
    cout<<value(sum)<<" "; }        // must be a Cartesian basis vector
  cout<<"\n";
}

