#ifndef VBF_mat_RR__H
#define VBF_mat_RR__H

#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ.h>
#include "vbf_vec_RR.h"

NTL_CLIENT

class vbf_mat_RR: public mat_RR {

  long IsIdent(const mat_RR& A, long n);

  int IsNotDefined(const mat_RR& a);

  void convolution(mat_RR& X, const mat_RR& A, const mat_RR& B);

  void Kronecker(mat_RR& X, const mat_RR& A, const mat_RR& B);

  RR maxvalue_abs(const mat_RR& X);
  // Calculate the maximum absolute value of the matrix X without
  // taking into account the first column and the first row

}; // end class vbf_mat_RR

void Kronecker(mat_RR& X, const mat_RR& A, const mat_RR& B)
{
   long n = A.NumRows();  
   long m = A.NumCols();  
   long p = B.NumRows();  
   long q = B.NumCols();  
   mat_RR C;
    
   X.SetDims(n*p, m*q);  
   C.SetDims(p,q);
  
   long i, j, k, l;  
   for (i = 0; i < n; i++)   
   {
      for (j = 0; j < m; j++)  
      {
         mul(C, A[i][j], B);  
         for (k = 0; k < p; k++)  
         {
            for (l = 0; l < q; l++)  
            {
               X[k+i*p][l+j*q] = C[k][l];
            }         
         }
      }
   }
}


void convolution(mat_RR& X, const mat_RR& A, const mat_RR& B)
{ 
   long n = A.NumRows();
   long m = A.NumCols();
   long q = B.NumCols();
 
   if (B.NumRows() != n)
      Error("matrix convolution: dimension mismatch");
 
   long cols = m*q;
   X.SetDims(n,cols);
   X[0][0] = n;
 
   long i, j, pos_a, pos_b;
   int bits_m, bits_q, bits_cols;
   mat_RR ATr, BTr, XTr;
   vec_GF2 x, a, b;
  
   ATr = transpose(A);
   BTr = transpose(B);
   XTr = transpose(X);
   bits_m = logtwo(m);
   bits_q = logtwo(q);
   a.SetLength(bits_m);
   b.SetLength(bits_q);
   bits_cols = bits_m+bits_q;

   for (i = 1; i < cols; i++)
   {
      x = to_vecGF2(i,bits_cols);
      for (j = 0; j < bits_m; j++)
      {
         a[j] = x[j];
      }
      for (j = 0; j < bits_q; j++)
      {
         b[j] = x[j+bits_m];
      }
      if (IsZero(a))
      {
         pos_b = conv_long(b);
         XTr[i] = BTr[pos_b];
      }
      else if (IsZero(b))
      {
         pos_a = conv_long(a);
         XTr[i] = ATr[pos_a];
      }
      else
      {
         pos_a = conv_long(a);
         pos_b = conv_long(b);
         convolution(XTr[i],ATr[pos_a],BTr[pos_b]);
      }
   }
   X = transpose(XTr);
} 

RR maxvalue_abs(const mat_RR& X)
{
   RR   max, temp;
   long i, j;
   long n = X.NumRows();
   long m = X.NumCols();

   max = 0;
   for (i = 1; i < n; i++)
   {
      for (j = 0; j < m; j++)
      {
         if (i != 0 || j != 0)
         {
             abs(temp,X[i][j]);
             if (temp > max) max = temp;
         }
      }
   }

   return max;
}

#endif
