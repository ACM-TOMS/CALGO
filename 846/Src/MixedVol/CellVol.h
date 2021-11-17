#ifndef _Cell_Vol_
#define _Cell_Vol_  

#include <iostream>

int gcd(int r1,int r2, int& k, int& l)
{
  k = 0,  l = 1;
  int xx = 1, yy = 0, q = r1/r2, r3 = r1%r2;
  while( r3 )
    {
      int tempxx = xx - q*k, tempyy = yy - q*l;
      xx = k; yy = l; k  = tempxx; l  = tempyy; r1 = r2; r2 = r3; q  = r1/r2; r3 = r1%r2;
    }
  return ( r2 < 0 )? (k = -k, l = -l, r2 = -r2) : r2;
}

void CellVol
   (
   int &nVar, int &nSpt, int** Spt, int* SptType,
   int* CurCell, int &MixedVol
   )

{
   int i, i1, j, j1, k, l, m, n, d;

   d = -1;
   n = -1;
   for ( i=0; i<nSpt; i++ )
   {
      l = CurCell[++d];
      for ( j=0; j<SptType[i]; j++ )
      {
         m = CurCell[++d];
         ++n;
         for ( k=0; k<nVar; k++ )  iwork[n][k] = Spt[m][k] - Spt[l][k];
      }
   }

   for ( i=0; i<nVar-1; i++ )
   {
      j = i;
      while ( j<nVar && iwork[j][i] == 0 )
         ++j;
      
      if ( j == nVar  )
         cerr << "Error code V001\n"; // volume of a cell = 0
      else
      {
         if( j != i )  
            std::swap(iwork[i],iwork[j]);
     
         for ( j=i+1; j<nVar; ++j)
            if( iwork[j][i] != 0 )
            {
               d = gcd(iwork[i][i], iwork[j][i], k, l);
               m = -iwork[j][i]/d;
               n = iwork[i][i]/d;
               for ( j1=i; j1<nVar; j1++)
               {
                  i1 = iwork[i][j1];
                  iwork[i][j1] = k*i1 + l*iwork[j][j1];
                  iwork[j][j1] = m*i1 + n*iwork[j][j1];
               }
            }
      }
   }
  
   d = iwork[nVar-1][nVar-1];
   for ( i=0; i<nVar-1; i++ )  d *= iwork[i][i];

   if (d < 0)
      MixedVol -= d;
   else
      MixedVol += d;

   return;
} 
 
#endif
