#ifndef VBF_mat_GF2__H
#define VBF_mat_GF2__H

#include <vector>
#include <NTL/mat_GF2.h>
#include <NTL/vec_long.h>
#include <algorithm>
#include "vbf_vec_GF2.h"
#include "vbf_ZZ.h"
#include "vbf_tools.h"

NTL_CLIENT

class vbf_mat_GF2: public mat_GF2 {

  int IsNotDefined(const mat_GF2& a);

  void concat_mat_GF2(mat_GF2& X, const mat_GF2& A, const mat_GF2& B);
  inline mat_GF2 concat_mat_GF2(const mat_GF2& A, const mat_GF2& B)
     { mat_GF2 X; concat_mat_GF2(X, A, B); NTL_OPT_RETURN(mat_GF2, X); }

  void directsum_mat_GF2(mat_GF2& X, const mat_GF2& A, const mat_GF2& B);

  void composition_mat_GF2(mat_GF2& X, const mat_GF2& A, const mat_GF2& B);
  // X = A comp B being A

  mat_GF2 rev(const mat_GF2& X, int n, int m);
  // Obtain Truth Table from ANF Table

  void juxtapose(mat_GF2& X, const mat_GF2& A, const mat_GF2& B);
  // X = A | B

  void reverse(mat_GF2& B, const mat_GF2& A);
  // B = reverse of A rows

  int IsConstant(const mat_GF2& A);

  mat_GF2 matGF2_seq(long rows, long cols);
  // Returns a matrix whose rows are the binary representation of 0,1,2,...rows-1

  void print(NTL_SNS ostream& s, const mat_GF2& a);

  void project(mat_GF2& X, const mat_GF2& A, const vec_long& b);

  int IsLinear(mat_GF2& X, const mat_GF2& A);
// If A is linear it returns 1 and the generator matrix X, else it returns 0.

  void mat_GF2tovec_long(vec_long& x, const mat_GF2& a);

  void conv_longtomat_GF2(mat_GF2& x, const vec_long& a);

  int aibf(const vec_GF2& f, int n, int d);

}; // end class vbf_mat_GF2

int IsNotDefined(const mat_GF2& a)
{
   long n = a.NumRows();

   if (n > 0) return 0;
  
   return 1;
}

void concat_mat_GF2(mat_GF2& X, const mat_GF2& A, const mat_GF2& B)
{
   long n = A.NumRows();  
   long m = A.NumCols();  
   long p = B.NumRows();  
   long q = B.NumCols(); 
   long row;   
   mat_GF2 C;
  
   X.SetDims(n*p, m+q);  
   C.SetDims(p,q);
  
   long i, j, k;  
  
   for (i = 1; i <= n; i++) {  
      for (j = 1; j <= p; j++) { 
	 row = j+(i-1)*p;
	 for (k = 1; k <= m; k++) { 
	    X(row,k)=A(i,k);
         }
	 for (k = 1; k <= q; k++) { 
	    X(row,k+m)=B(j,k);
         }
      }
   }  
}

void directsum_mat_GF2(mat_GF2& X, const mat_GF2& A, const mat_GF2& B)
{
   long n = A.NumRows();  
   long m = A.NumCols();  
   long p = B.NumRows();  
   long i,j,k,row;  

   if (m != B.NumCols()) Error("directsum: Image dimension mismatch");
  
   X.SetDims(n*p, m);  

   int numbits = logtwo(n);
  
   for (i = 0; i < n; i++) {  
      for (j = 0; j < p; j++) {  
	 for (k = 0; k < m; k++) {
	    row = (i<<numbits)+j;
	    X[row][k] = A[i][k]+B[j][k];
	 }
      }
   }  
}

// A o B
void composition_mat_GF2(mat_GF2& X, const mat_GF2& A, const mat_GF2& B)
{
   long ra = A.NumRows();  
   long ca = A.NumCols();  
   long rb = B.NumRows();  
   long cb = B.NumCols();
   long r, i, pos; 
   vec_GF2 v;

   r = logtwo(ra);
   if (cb != r) Error("composition: dimensions mismatch");
  
   X.SetDims(rb, ca); 

   for (i = 0; i < rb; i++) 
   {
      v = B[i];
      pos = conv_long(v);
      X[i] = A[pos];
   } 

}

mat_GF2 rev(const mat_GF2& X, int n, int m)
{
   mat_GF2 Y, Z;
   long    k, mi, numrows;
   int     i, j;

   numrows = (1 << n);

   Y.SetDims(numrows, m);
   Z.SetDims(numrows, m);
   Y = X;
  
   for (j = 0; j < m; j++) 
   {
      mi = 1;
      for (i = 0; i < n; i++) 
      {
         for (k = 0; k < numrows; k++)
           if ((k >> i) % 2) Z[k][j] = Y[k-mi][j] + Y[k][j];
           else Z[k][j] = Y[k][j];
         for (k = 0; k < numrows; k++) Y[k][j]= Z[k][j];
         mi *= 2;
      }   
   }
  
   return Y;
}

void juxtapose(mat_GF2& X, const mat_GF2& A, const mat_GF2& B)
{
    long i, j;  
    int an = A.NumRows();
    int am = A.NumCols();
    int bn = B.NumRows();  
    int bm = B.NumCols();

    if (an != bn)   
       Error("juxtapose: Row dimension mismatch");  
         
    X.SetDims(an,(am+bm));
    for (i = 0; i < an; i++)  
    {
       for (j = 0; j < am; j++)
       {
   	  X[i][j]=A[i][j];
       }	
       for (j = am; j < (am+bm); j++)
       {
   	  X[i][j]=B[i][j-am];
       }	
    }  

}

int IsConstant(const mat_GF2& A)
{	
   long i, j;
     
   for (i = 0; i < A.NumRows(); i++) 
   { 
      for (j = 0; j < (A.NumCols()-1); j++)
      {
      	 if (A[i][j] != A[i][j+1]) return 0;
      }
   }
   return 1;
}

void reverse(mat_GF2& B, const mat_GF2& A)
// B = reverse of A rows
{
   long i,n,m;

   n = A.NumRows();
   m = A.NumCols();
   B.SetDims(n,m);

   for (i= 0;i < n;i++) 
   {
      reverse(B[i],A[i]);
   }

}

mat_GF2 matGF2_seq(long rows, long cols)
{
   mat_GF2 A;
   long i;

   A.SetDims(rows,cols);
   for (i = 0; i < rows; i++)
   {
      A[i] = to_vecGF2(i,cols);
   }
   return A;
}
	
void print(NTL_SNS ostream& s, const mat_GF2& a)  
{  
   long i, j;  
      
   for (i = 0; i < a.NumRows(); i++) 
   {
      for (j = 0; j < a.NumCols(); j++)
      {
         s << " " << a[i][j];
      }   
      s << "\n";	 
   }
}  

void project(mat_GF2& X, const mat_GF2& A, const vec_long& b)
{
    long i;  
    int an = A.NumRows();
    int am = A.NumCols();
    int bn = b.length();  
    mat_GF2 ATr, XTr;

    if (am < bn)   
       Error("project: Excesive coordinates to be projected");  
         
    X.SetDims(an,bn);
    ATr = transpose(A);
    XTr = transpose(X);

    for (i = 0; i < bn; i++)  
    {
       XTr[i] = ATr[b[i]];
    }  
    X = transpose(XTr);
  
}

// If A is linear it returns 1 and the generator matrix X, else it returns 0.
int IsLinear(mat_GF2& X, const mat_GF2& A)
{
    long i;
    int spacen = A.NumRows();
    int m = A.NumCols();
    int n;
    vec_GF2 bin;
    mat_GF2 B;

    n = logtwo(spacen);
    X.SetDims(n,m);
    B.SetDims(spacen,m);

    X[0]=A[1];

    for (i = 1; i < n; i++)
	X[i] = A[1 << i];

    for (i = 0; i < spacen; i++)
    {
        bin = to_vecGF2(i,n);
        mul(B[i],bin,X);
    }
   
    if (A == B)
    {
	return 1;
    } else {
	return 0;
    } 
}

void mat_GF2tovec_long(vec_long& x, const mat_GF2& a)
{
   long n = a.NumRows();

   if (n == 0) {
      x.SetLength(0);
      return;
   }

   long i;

   x.SetLength(n);
   for (i = 0; i < n; i++)
      x[i] = conv_long(a[i]);
} 

void conv_longtomat_GF2(mat_GF2& x, const vec_long& a)
{
   long n = a.length();

   if (n == 0) {
      x.SetDims(0,0);
      return;
   }

   long i, max = 0;
   int m = 1;

   for (i = 0; i < n; i++)
      if (a[i] > max) max = a[i];

   while (conv_long(to_vecGF2(max,m)) != max)
      m++;

   x.SetDims(n,m);
   for (i = 0; i < n; i++)
      x[i] = to_vecGF2(a[i],m);
}

void M(mat_GF2& X, const int& n, const int& deg, const vec_GF2& f, const long& t)
{
   int i,j;
   long x,rows,spacen,num,cont=1;
   long *v = NULL;
   NTL::ZZ a=to_ZZ(1), b, c;
   NTL::mat_GF2 T;
   vector<long> support;

   spacen = 1 << n;
   NTL::mat_GF2 A(INIT_SIZE, spacen,1);

   rows = weight(f);
   for (i = 0; i < f.length(); i++)
   {
      if (f[i] == 1) support.push_back(i);
   }
 
   X.SetDims(rows,t);

   for (i = 0; i < rows; i++)
   {
      X[i][0] = 1;
   }
  
   for (i = 1; i <= deg; i++)
   {
      num = to_long(numofweight(n,i));
      v = (long *) malloc(num * sizeof(long));
      vectors_weight(v, n, i);

      vector<long> vsorted (v,v+num);
      vector<long>::iterator it;
      std::sort(vsorted.begin(),vsorted.end());

      for (it=vsorted.begin(); it!=vsorted.end(); ++it)
      {
         clear(A);
         x = *it;
         A[x][0] = 1;
         T = rev(A, n, 1);

         for (j = 0; j < rows; j++)
         {
	    X[j][cont] = T[support[j]][0];
	 }
	 cont++;
      }
   }
    
}

int aibf(const vec_GF2& f, int n, int d)
{
   int  deg;
   long r1,r2,r,t;
   vec_GF2 fn;
   mat_GF2 X;

   t = 1;

   for (deg = 1; deg <= d-1; deg++)
   {
       t = t + Combination(n,deg);

       M(X,n,deg,f,t);
       r1 = gauss(X);
       NTL::negate(fn,f);
       M(X,n,deg,fn,t);
       r2 = gauss(X);

       if (r1 < r2)
       {
          r = r1;
       } else {
          r = r2;
       }

       if (r < t) return deg;
   }

   return d;
}


#endif
