#ifndef _Lower_Triangular_
#define _Lower_Triangular_

void LowerTriangular
   (
   double** A, int jStrt, int jEnd, int na, 
   int &rnk, int* ib
   )

/*
   On exit
   rnk : the rank of that portion of the matrix.
   ib : the row indices of lower-triagular part

 * Purpose
 * =======
 
   Lower-Triangulize the rows from jStrt to jEnd of matirx A with na columns.
   rnk is the rank of that portion of A on exit.
   No row interchange is allowed.
*/

{
   int i, j, k, itmp;
   double dtmp;

   // Since the jStrt-th row of A is (1,0,...,0)
   ib[0] = jStrt;
   for (i=1; i<na; i++) ib[i] = -1;
   rnk = 1;

   k = jStrt+1;
   while ( rnk < na && k <= jEnd )
   {
      // to find the biggest component from A[k][rnk] to A[k][na-1]
      dtmp = 1.0e-13;
      itmp = -1;
      for ( i = rnk; i<na; i++ )
         if ( fabs(A[k][i]) > dtmp )
         {
            dtmp = fabs(A[k][i]);
            itmp =i;
         }

      if ( itmp >= 0)  // non-zero row
      {
         // To change aa(k,*) into e_itmp
         for (i=0; i<itmp; i++) A[k][i] /= A[k][itmp];
         for (i=itmp+1; i<na; i++) A[k][i] /= A[k][itmp];
         for (j=k+1; j<=jEnd; j++)
         {
            for (i=0; i<itmp; i++) A[j][i] -= A[k][i]*A[j][itmp];
            for (i=itmp+1; i<na; i++) A[j][i] -= A[k][i]*A[j][itmp];
            A[j][itmp] /= A[k][itmp];
	 }
         if ( itmp != rnk ) // interchange the rnk-th and itmp-th columns of A
         {
            for ( j=jStrt; j<=jEnd; j++ )
            {
               dtmp = A[j][rnk];
               A[j][rnk] = A[j][itmp];
               A[j][itmp] = dtmp;
            }
         }

         for ( i=0; i<rnk; i++ ) A[k][i] = 0.0e0;
	 A[k][rnk] = 1.0e0;
         for ( i=rnk+1; i<na; i++ ) A[k][i] = 0.0e0;

         ++rnk;
	 ib[rnk-1]=k;
      }
      ++k;
   }
}

#endif
