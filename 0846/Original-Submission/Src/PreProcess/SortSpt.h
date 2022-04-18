#ifndef _Sort_Supports_
#define _Sort_Supports_

void SortSpt
   (
      int &nVar, int &nSpt, int** Spt, int* Spt1stIdx
   )
/*
   nVar : the dimension of the supports.
   nSpt : the number of supports.
   Spt  : the supports; each row is a point.
   Spt1stIdx : the index of the 1st point of each support in Spt.
               Spt1stIdx[nSpt]-1 is the index of the last point of the supports.

   Purpose:
   ========

   To sort the points in each support of the nSpt supports 
   in decending lexicographic order.
*/

{
   int iSpt, j, k, idx;

   for (iSpt=0; iSpt<nSpt; iSpt++)
   {
      for (int i=Spt1stIdx[iSpt]; i<Spt1stIdx[iSpt+1]-1; i++)
      {
         idx = i;
         for (j=i+1; j<Spt1stIdx[iSpt+1]; j++)
	 {
            for (int k=0; k<nVar; k++)
	    {
               if (Spt[j][k] < Spt[idx][k])
	       {
	          break;  // LexiGt = false
	       }
               else
	       {
                  if (Spt[j][k] > Spt[idx][k])
		  {
                     idx = j;
		     break; // LexiGt = true;
		  }
	       }
	    } 
	 }
         if (idx>i)
	 {
            for (j=0; j<nVar; j++)
	    {
               k = Spt[i][j];
               Spt[i][j] = Spt[idx][j];
               Spt[idx][j] = k;
	    }
	 }
      }
   }
}

#endif
