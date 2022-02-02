#ifndef VBF_mat_ZZ__H
#define VBF_mat_ZZ__H

#include <NTL/mat_ZZ.h>
#include <NTL/mat_GF2.h>
#include <NTL/RR.h>
#include "vbf_vec_ZZ.h"

NTL_CLIENT

class vbf_mat_ZZ: public mat_ZZ {

  int IsNotDefined(const mat_ZZ& a);

  void div(mat_ZZ& X, const mat_ZZ& A, const long& b_in); 
  // Divide all the elements of the matrix A by the number b_in

  inline mat_ZZ div(const mat_ZZ& A, const long& b_in)	
  { mat_ZZ X; div(X, A, b_in); return X; }

  ZZ maxvalue_abs(const mat_ZZ& X);
  // Calculate the maximum absolute value of the matrix X without
  // taking into account the first column

  ZZ maxvalue(const mat_ZZ& X);
  // Calculate the maximum value of the matrix X without
  // taking into account the first column

  mat_ZZ wt(const mat_ZZ& X, int n, int m);
  // Calculate Walsh Transform

  mat_ZZ invwt(const mat_ZZ& X, int n, int m);
  // Calculate Inverse Walsh Transform

  mat_ZZ lat(const mat_ZZ& X, int spacen, int spacem);
  // Calculate Linear Approximation table

  int WalshIsZero(const mat_ZZ& X, int n, int m, int t);
  // Find-out if the nxm Matrix X has the rows corresponding
  // to the vectors with weight t equal to zero

  int IsConstant(const mat_ZZ& X, int n, int m, int w);
  // Find-out if the nxm Matrix X has the rows corresponding
  // to the vectors with weight t as constant vectors

  mat_ZZ charfunct(const mat_GF2& T, int n, int m);
  // Calculate Characteristic function from the Truth table

  mat_ZZ sequence(const mat_GF2& T);
  // Calculate the sequence of the S-box with the Truth table T

  mat_GF2 truthtable(const mat_ZZ& C, int n, int m);
  // Calculate the Truth Table from the Characteristic function

  void convolution(mat_ZZ& X, const mat_ZZ& A, const mat_ZZ& B);

  void print(NTL_SNS ostream& s, const mat_ZZ& a);

  mat_GF2 to_tt(const mat_ZZ& S);
  // Calculate the Truth Table from the Sbox representation

  void Kronecker(mat_ZZ& X, const mat_ZZ& A, const mat_ZZ& B);
  // Calculate the Kronecker product of A and B 

  mat_ZZ Sylvester(int n);
  // Calculate the Sylvester-Hadamard matrix of order 2^n

  mat_GF2 seq_to_tt(const mat_ZZ& S);
  // Calculate the sequence from the truth table

  void energy(mat_ZZ& X, const mat_ZZ& A); 
  // X(i,j) = A(i,j)^2

  inline mat_ZZ energy(const mat_ZZ& A)	
  { mat_ZZ X; energy(X, A); return X; }

  void num_zeros(ZZ& c, const mat_ZZ& A);

  inline ZZ num_zeros(const mat_ZZ& A)	
  { ZZ c; num_zeros(c, A); return c; }

  void convol_column(mat_ZZ& X, const mat_ZZ& A, const mat_ZZ& B);

  void project(mat_ZZ& X, const mat_ZZ& A, const vec_long& b);

  ZZ maxall(const mat_ZZ& X);

  ZZ maxall_abs(const mat_ZZ& X);

}; // end class vbf_mat_ZZ

int IsNotDefined(const mat_ZZ& a)
{
   long n = a.NumRows();

   if (n > 0) return 0;
  
   return 1;
}

void div(mat_ZZ& X, const mat_ZZ& A, const long& b_in)  
{
   long b = b_in;
   long n = A.NumRows();
   long m = A.NumCols();

   X.SetDims(n, m);

   long i, j;	 
   for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
         X[i][j]= A[i][j]/b;
} 

ZZ maxvalue_abs(const mat_ZZ& X)
{
   ZZ   max, temp;
   long i, j;
   long n = X.NumRows();
   long m = X.NumCols();

   max = 0;
   for (i = 0; i < n; i++)
   {
      for (j = 1; j < m; j++) 
      {
         abs(temp,X[i][j]);
         if (temp > max) max = temp;
      }	
   }  

   return max;
}

ZZ maxvalue(const mat_ZZ& X)
{
   ZZ   max;
   long i, j;
   long n = X.NumRows();
   long m = X.NumCols();

   max = 0;
   for (i = 0; i < n; i++)
   {
      for (j = 0; j < m; j++) 
      {
         if (i != 0 && j != 0)
         {
            if (X[i][j] > max) max = X[i][j];
	 }
      }	
   }  

   return max;
}

ZZ maxall(const mat_ZZ& X)
{
   ZZ   max;
   long i, j;
   long n = X.NumRows();
   long m = X.NumCols();

   max = 0;
   for (i = 0; i < n; i++)
   {
      for (j = 0; j < m; j++)
      {
         if (X[i][j] > max) max = X[i][j];
      }
   }

   return max;
}

ZZ maxall_abs(const mat_ZZ& X)
{
   ZZ   max,temp;
   long i, j;
   long n = X.NumRows();
   long m = X.NumCols();

   max = 0;
   for (i = 0; i < n; i++)
   {
      for (j = 0; j < m; j++)
      {
         abs(temp,X[i][j]);
         if (temp > max) max = temp;
      }
   }

   return max;
}

mat_ZZ wt(const mat_ZZ& X, int n, int m)
{
   mat_ZZ 	Y;
   long i, j, k, l, mi, mia, a, b, numrows, numcolumns, numtotal;

   numrows = (1 << n);
   numcolumns = (1 << m);
   numtotal = (1 << (n+m));
   
   mat_ZZ Z(INIT_SIZE, numrows, numcolumns);
  
   Y = X;
   
   mi = 1;
   for (i = 0; i < (n+m); i++) 
   {
      for (j = 0; j < numtotal; j++)
      {
      	k = j / numcolumns;
      	l = j % numcolumns;
        if ((j >> i) % 2)
        {
           if (mi <= l)
           { 
              a = k;
              b = l-mi;
           }
           else 
           {
              mia = mi-l;
              b = mia % numcolumns;
              if (b == 0)
              {
              	 a = k-(mia / numcolumns);
              }
              else
              {
              	 a = k-(mia / numcolumns)-1;
              	 b = numcolumns-b;
              }	                   
           }   
           Z[k][l] = Y[a][b] - Y[k][l];
        }  
        else
        {
           if ((l+mi) < numcolumns)
           { 
              a = k;
              b = l+mi;
           }
           else 
           {
              mia = mi-(numcolumns-l);
              b = mia % numcolumns;
              a = k + (mia / numcolumns) +1;                 
           }   
           Z[k][l] = Y[k][l] + Y[a][b]; 
        }       	
      }	
      Y = Z;
      mi *= 2;
   }

   return Y;
}  

mat_ZZ invwt(const mat_ZZ& X, int n, int m)
{
   mat_ZZ	Winv, C;
   long 	numrows, numcolumns, d;
         
   Winv = wt(X, n, m);
   numrows = (1 << n);
   numcolumns = (1 << m);

   d = numrows * numcolumns;
   div(C, Winv, d);

   return C;
}  

mat_ZZ lat(const mat_ZZ& X, int spacen, int spacem)
{
   mat_ZZ  LAT;
   long	   i, j;
   
   LAT.SetDims(spacen,spacem);	

   for (i = 0; i < spacen; i++)
   {
      for (j = 0; j < spacem; j++) 
      {
         LAT[i][j] = X[i][j]*X[i][j];
      }	
   }  
   return LAT;	 	    
}

int WalshIsZero(const mat_ZZ& X, int n, int m, int t)
{
   long 	num, i, spacen, spacem;
   long 	*v = NULL;
   vec_ZZ	F;
   
   spacen = (1 << n);
   spacem = (1 << m);
   
   if (t == 0) 
   {
      F.SetLength(spacem);
      F = X[0];
      if (F[0] == spacen)
      {
        F[0] -= F[0];
        if (!IsZero(F)) return 0;
      }
   } else
   {  
      num = to_long(numofweight(n, t));
      v = (long *) malloc(num * sizeof(long));
      vectors_weight(v, n, t);

      for (i = 0; i < num; i++) 
      { 
        if (!IsZero(X[v[i]])) return 0;
      }
   }         
   return 1;
}

int IsConstant(const mat_ZZ& X, int n, int m, int w)
{
   long 	num, i, j, spacem;
   long 	*v = NULL;
   
   num = to_long(numofweight(n, w));
   v = (long *) malloc(num * sizeof(long));
   vectors_weight(v, n, w);
   spacem = (1 << m);
   
   for (i = 0; i < num; i++) 
   { 
      for (j = 0; j < spacem-1; j++)
      {
      	 if (X[v[i]][j] != X[v[i]][j+1]) return 0;
      }
   }
   return 1;
}

mat_ZZ charfunct(const mat_GF2& T, int n, int m)
{
   long i, j, y, numrows, numcolumns;

   numrows = (1 << n);
   numcolumns = (1 << m);
   
   mat_ZZ X(INIT_SIZE, numrows, numcolumns); 

   for (i = 0; i < numrows; i++) 
   {                 
      y = 0;                        
      for (j = 0; j < m; j++) 
      {
    	 if (T[i][m-1-j] == 1) 
    	 {
    	    y = y + (1 << j);
    	 }   	
      }
      X[i][y] = 1;              
   }

   return X;	
}	

mat_ZZ sequence(const mat_GF2& T)
{
   long i, j, numrows = T.NumRows(), numcols = T.NumCols();
   
   mat_ZZ X(INIT_SIZE, numrows, numcols); 

   for (i = 0; i < numrows; i++) 
   {                 
      for (j = 0; j < numcols; j++) 
      {
    	 if (T[i][j] == 1) 
    	 {
    	    X[i][j] = -1;
    	 } else {
    	    X[i][j] = 1; 
         }  	
      }
                    
   }

   return X;	
}	

mat_GF2 seq_to_tt(const mat_ZZ& S)
{
   long i, j, numrows = S.NumRows(), numcols = S.NumCols();
   
   mat_GF2 X(INIT_SIZE, numrows, numcols); 

   for (i = 0; i < numrows; i++) 
   {                 
      for (j = 0; j < numcols; j++) 
      {
    	 if (S[i][j] == 1) 
    	 {
    	    X[i][j] = 0;
    	 } else {
    	    X[i][j] = 1; 
         }  	
      }
                    
   }

   return X;	
}	

mat_GF2 truthtable(const mat_ZZ& C, int n, int m)
{
   long    i, j, numrows;
   mat_GF2 T;
   vec_GF2 v;   

   numrows = (1 << n);
      
   T.SetDims(numrows,m);
   v.SetLength(m);
   
   // Calculate the truthtable
   for (i = 0; i < numrows; i++) 
   {         
//cout << C[i] << endl;
      if (IsVecCharFunct(C[i],j))
      {
         v = to_vecGF2(j, m);
         T[i] = v;
      }
      else
      {
//       Error("truthtable: Characteristic Function not valid");
         T.kill();
	 return T;
      } 	    
   }
   return T;	
}	

void convolution(mat_ZZ& X, const mat_ZZ& A, const mat_ZZ& B)  
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
   mat_ZZ ATr, BTr, XTr;
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

void print(NTL_SNS ostream& s, const mat_ZZ& a)  
{  
   long spacen = a.NumRows();  
   long m = a.NumCols();
   long i, j, k, num;  
   int n;
   long *v = NULL;
      
   n = logtwo(spacen);

   s << "0";
   for (j = 1; j < m; j++)
   {
      s << " " << to_RR(a[0][j])/to_RR(a[0][0]);
   }     
   s << "\n";
   for (i = 1; i < spacen-1; i++) 
   {
      num = to_long(numofweight(n, i));
      v = (long *) malloc(num * sizeof(long));
      vectors_weight(v, n, i);
         
      for (j = 0; j < num; j++)
      {
         s << v[j];
      	 for (k = 1; k < m; k++)
      	 {
            s << " " << to_RR(a[v[j]][k])/to_RR(a[0][0]);
         }   
         s << "\n";
      }     
   }
   s << (spacen-1);
   for (j = 1; j < m; j++)
   {
      s << " " << to_RR(a[spacen-1][j])/to_RR(a[0][0]);
   }     
   s << "\n";
}  

mat_GF2 to_tt(const mat_ZZ& S)
{
   long i, j, k, r, c, row, cont, impar, pos;
   long numcols = S.NumCols();
   long numrows = S.NumRows();
   mat_GF2 T;
   vec_GF2 vr, vc, vt, vs;

   r = logtwo(numrows);
   c = logtwo(numcols);       
   T.SetDims(numrows*numcols,c);
   vs.SetLength(r+c);

   // Calculate the truthtable
   for (i = 0; i < numrows; i++) 
   {         
      for (j = 0; j < numcols; j++) 
      {         
         vr = to_vecGF2(i,r);
         vc = to_vecGF2(j,c);
         vt = to_vecGF2(to_long(S[i][j]),c);
         if (r > 0) {
            impar = 1;
            cont = 0;
            while (cont < r) {
               pos = (cont/2)+(cont%2);
               if (impar == 1) { 
                  vs[pos] = vr[cont];
                  impar = 0;
               } else {   
                  vs[r+c-pos] = vr[cont];
                  impar = 1;
               }
               cont++;
            }
            if (impar == 0) {
               pos = pos+1;
            } else {
	       pos = (cont/2)+(cont%2);
            }
            for (k = 0; k < c; k++) 
            {         
               vs[k+pos] = vc[k];
            }
         }
	 else {
            vs = vc;
         }
         row = conv_long(vs);
         T[row] = vt;
      }
   }

   return T;	
}	

void Kronecker(mat_ZZ& X, const mat_ZZ& A, const mat_ZZ& B)
{
   long n = A.NumRows();  
   long m = A.NumCols();  
   long p = B.NumRows();  
   long q = B.NumCols();  
   mat_ZZ C;
    
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

mat_ZZ Sylvester(int n)
{
   mat_ZZ A, B, X;
   int i;
   long r;

   A.SetDims(1,1);
   A[0][0] = 1;
   B.SetDims(2,2);
   B[0][0] = 1;
   B[0][1] = 1;
   B[1][0] = 1;
   B[1][1] = -1;

   for (i = 0; i < n; i++)   
   {
      Kronecker(X,B,A);
      r = X.NumRows();
      A.SetDims(r,r);
      A = X;     
   }   
   
   return X;
}

void energy(mat_ZZ& X, const mat_ZZ& A)
{
   long n = A.NumRows();  
   long m = A.NumCols();  
    
   X.SetDims(n, m);  
  
   long i, j;  
   for (i = 0; i < n; i++)   
   {
      for (j = 0; j < m; j++)  
      {
         X[i][j] = A[i][j]*A[i][j];  
      }
   }

}

void num_zeros(ZZ& c, const mat_ZZ& A)
{
   long n = A.NumRows();  
   long m = A.NumCols();  
   
   c = 0; 
   long i, j;  
   for (i = 0; i < n; i++)   
   {
      for (j = 0; j < m; j++)  
      {
         if (A[i][j] == 0) c++;  
      }
   }

}

void convol_column(mat_ZZ& X, const mat_ZZ& A, const mat_ZZ& B)  
{  
   long n = A.NumRows();  
   long m = A.NumCols();
  
   if (B.NumRows() != n || B.NumCols() != m)   
      Error("convolution of columns: dimension mismatch");  
  
   X.SetDims(n, m);  
  
   long i; 
   mat_ZZ ATr, BTr, XTr;
   
   ATr = transpose(A);
   BTr = transpose(B);
   XTr = transpose(X);

   for (i = 0; i < m; i++)   
   {
      convolution(XTr[i],ATr[i],BTr[i]);
   }
   X = transpose(XTr);
}  

void project(mat_ZZ& X, const mat_ZZ& A, const vec_long& b)  
{  
   long ra = A.NumRows();  
   long ca = A.NumCols();
   long i, j, pos_a, am, cx, xm = b.length();
   mat_ZZ ATr, XTr;
   vec_GF2 x, a;

   am = logtwo(ca);
   if (am < xm)   
       Error("project: Excesive coordinates to be projected");  
         
   cx = (1 << xm);
   X.SetDims(ra,cx);
   ATr = transpose(A);
   XTr = transpose(X);
   a.SetLength(am);
   
   for (i = 0; i < cx; i++)   
   {
      x = to_vecGF2(i,xm);
      for (j = 0; j < xm; j++)   
      {
         a[b[j]] = x[j];
      }
      pos_a = conv_long(a);
      XTr[i] = ATr[pos_a];
      clear(a); 
   }
   X = transpose(XTr);
}  

#endif
