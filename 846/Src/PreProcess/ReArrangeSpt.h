#ifndef _Re_Arrange_Support_
#define _Re_Arrange_Support_

void ReArrangeSpt
   (
      int &nVar, int** Spt, int* Spt1stIdx, int** SptVtx, int* SptVtx1stIdx,
      int* SptType, int*NuIdx2OldIdx
   )
/*
 * nVar : dimension of the support
 * Spt : the support
 * Spt1stIdx : the position of the first point of each support in Spt.
 * SptVtx : the smaller support after deleting all non-vertex points.
 * SptVtx1stIdx : the position of the first point of each support in SptVtx.
 * SptType : On input: encoded with type information, 
 *                     SptType[i]=-1 means the i-th Spt is the base support;
 *                     SptType[i]=j means the i-th Spt = j-th Spt.
 *           On exit: the type of the support.
 * NuIdx2OldIdx : positions of the points in SptVtx to positions in Spt
*/

/* Purpose:
   =======

   To re-arrange the supports Spt and SptVtx according to the type determined.
  
   For example: A1, A2, A3 are same and A4 and A5 are same, then
   we re-arrange the supports into A1, A4, A1, A1, A4.
*/

{
   int i, i1, icnt, iSpt, ith, j;

   int* iTmp = new int [nVar];

   i = 0;
   iSpt = -1;
   while (i<nVar)
   {
      if (SptType[i] < 0)
      {
         ++iSpt;
         iTmp[iSpt] = i;
      }
      ++i;
   }
   i = 0;
   ith = -1; 
   while (i<nVar)
   {
      if (SptType[i] < 0)
      {
	 ++ith;
         icnt = 1;
         for (j=i+1; j<nVar; j++)
         {
            if (SptType[j] == i)
	    {
               ++iSpt;
               ++icnt;
               iTmp[iSpt] = j;
	    }
	 }
         SptType[ith] = icnt;
      }
      i = i + 1;
   }
   // SptType[i] = k : # of supports which are the same as the i-th base support.

   // To rearrange Spt, Spt1stIdx

   int** SptCopy = new int* [Spt1stIdx[nVar]];
   SptCopy[0] = new int [nVar*Spt1stIdx[nVar]];
   for (i=1; i<Spt1stIdx[nVar]; i++) SptCopy[i] = SptCopy[i-1] + nVar;
   int* Spt1stIdxCopy = new int [nVar+1];

   for (j=0; j<Spt1stIdx[nVar]; j++)
   {
      for (i=0; i<nVar; i++) SptCopy[j][i] = Spt[j][i];
   }

   icnt = -1;
   for (i=0; i<nVar; i++)
   {
      for (j=Spt1stIdx[iTmp[i]]; j<Spt1stIdx[iTmp[i]+1]; j++)
      {
         ++icnt;
         for (i1=0; i1<nVar; i1++) Spt[icnt][i1] = SptCopy[j][i1];
      }
      Spt1stIdxCopy[i+1] = icnt+1;
   } 
   Spt1stIdxCopy[0] = 0;

   // To rearrange NuIdx2OldIdx

   int* Nu2OldCopy = new int [Spt1stIdx[nVar]];

   for (j=0; j<SptVtx1stIdx[nVar]; j++) Nu2OldCopy[j] = NuIdx2OldIdx[j];

   icnt = -1;
   for (i=0; i<nVar; i++)
   {
      for (j=SptVtx1stIdx[iTmp[i]]; j<SptVtx1stIdx[iTmp[i]+1]; j++)
      {
         ++icnt;
         NuIdx2OldIdx[icnt] = Nu2OldCopy[j] - (Spt1stIdx[iTmp[i]] - Spt1stIdxCopy[i]);
      }
   }
 
   for (i=0; i<=nVar; i++) Spt1stIdx[i] = Spt1stIdxCopy[i];
 
   // To rearrange SptVtx, SptVtx1stIdx
 
   for (j=0; j<SptVtx1stIdx[nVar]; j++)
      for (i=0; i<nVar; i++) SptCopy[j][i] = SptVtx[j][i];

   icnt = -1;
   for (i=0; i<nVar; i++)
   {
      for (j=SptVtx1stIdx[iTmp[i]]; j<SptVtx1stIdx[iTmp[i]+1]; j++)
      {
         ++icnt;
         for (i1=0; i1<nVar; i1++) SptVtx[icnt][i1] = SptCopy[j][i1];
      }
      Spt1stIdxCopy[i+1] = icnt+1;
   }
   SptVtx1stIdx[0] = 0;
   for (i=1; i<=nVar; i++) SptVtx1stIdx[i] = Spt1stIdxCopy[i];

   delete [] iTmp;
   delete [] Nu2OldCopy;
   delete [] Spt1stIdxCopy;
   delete [] SptCopy[0];
   delete [] SptCopy;

   return;
}

#endif
