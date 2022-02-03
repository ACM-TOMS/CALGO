#ifndef _Revised_Simplex_Method_
#define _Revised_Simplex_Method_

#include <fstream>
#include <iostream>

void RSimplex
   (
    double** A, double* b, int jStrt, int jEnd, int m,
    double* xb, int* ib, double** binv, int &info
   )
/*
   A : Active part is row jStrt to row jEnd row  and  column 1 to column m.
   b : Active part is elemen jStrt to Element jEnd
   xb : the initial feasible solution.
   binv : the inverse of the initial base matrix.
   ib : the indices of the basics.
   info : position of the variable epsilon when info>=0.

   Purpose
   =======

   This is the revised simplex method for solving
                 min e
   under constraints
		 - e      <=0
                 - e + By <= b
   where (e,y) are the variables of the LP. -e <=0 corresponds to
   the jStrt-th row of A. The active part of A is from the jStrt-th
   row to jEnd-th row.
*/

{
   int i, j, k, ell;
   bool ibrk;
   double sigb, sum, vkmin; 
   double eps = 1.0e-6;
   double tol = 1.0e-10;

   extern double* sb;   //double* sb = new double [m];

   info = -1;
   for (i=0; i<m; i++)
   {
      if (ib[i] == jStrt)
      {
         info = i;
         break;
      }
   }

   while(1)
   {
      // To compute u and find k
      if (info == -1) // $e_{n+1}$ is not in the basis
      {
	 k = jStrt;
         vkmin = -1.0e0;
      }
      else  // $e_{n+1}$ is in the basis
      {
	 // To find vkmin=v_k=min{c_i+a_i^Tu | i\notin I_b}
         k = -1;
         vkmin = DBL_MAX;
         for (i=jStrt; i<=jEnd; i++)
	 {
            // To check if i is in ib
            ibrk = 0;
            for (j=0; j<m; j++)
            {
	       if (i == ib[j])
	       {
		  ibrk = 1;
		  break;
	       }
	    }
	    if (ibrk) continue;

            // If i is not in ib
            sum = 0.0e0;
            for (j=0; j<m; j++) sum += A[i][j]*binv[info][j];
            if (sum < vkmin)
	    {
	       vkmin = sum;
               k = i;
            }
	 }
         if (vkmin > -eps) return; // Have found the optimal solution
      }

      // k-th column will get into the basis
      // To form sb
      for (i=0; i<m; i++)
      {
         sum = 0.0e0;
         for (j=0; j<m; j++) sum += binv[i][j]*A[k][j];
         sb[i] = sum;
      }

      ell = 0;
      sigb = DBL_MAX;
      for (i=0; i<m; i++)
      {
	 if (sb[i] > eps && (xb[i]/sb[i]) < sigb) 
	 {
	    sigb = (xb[i]/sb[i]);
            ell = i;
         }
      }
      // ell-th row gets out of the basis and k-th row gets into the basis.
      // To find the new B^{-1}
      sum = 0.0e0;
      for (i=0; i<m; i++) sum += binv[ell][i]*A[k][i];
 
      if (fabs(sum) <= eps)
      {
        // cerr << "B^{-1} is singular" << endl;
         abort();
      }

      sum = 1.0e0/sum;
      for (i=0; i<m; i++) binv[ell][i] *= sum;

      for (j=0; j<m; j++)
      {
         if (j != ell)
	 {
            sum = 0.0e0;
            for (i=0; i<m; i++) sum += binv[j][i]*A[k][i];

            for (i=0; i<m; i++) binv[j][i] -= sum*binv[ell][i];
         }
      }

      // To form the new basic feasible solution
      for (i=0; i<m; i++)
      {
	 sum = 0.0e0;
         for (int j=0; j<m; j++) sum += binv[i][j]*b[j];
         xb[i] = sum;
      }

      if (ib[ell] == jStrt) info = 0;
      // the jStrt-th column is out of basis

      ib[ell] = k;
      if (k == jStrt)
      {
         info = ell;
         if (xb[info] > tol) return;  // epsilon is already positive.
      }
   }
}

#endif

