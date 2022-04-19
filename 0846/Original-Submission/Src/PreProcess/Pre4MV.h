#ifndef _Preprocess_for_Mixed_Volume_Supports_
#define _Preprocess_for_Mixed_Volume_Supports_

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <numeric>
#include <fstream>
#include <iostream>

#include "SortSpt.h"
#include "NonVertex.h"
#include "ReArrangeSpt.h"

void Pre4MV
   (
      int &nVar, int &nSpt, int* SptType,
      int** Spt, int* Spt1stIdx,
      int** SptVtx, int* SptVtx1stIdx,
      int* NuIdx2OldIdx
   )
/*
 * Purpose:
 * To preprocess the supports for mixed volume computation
*/

{
   int i, j, k, l, icnt, itmp;

   // To sort the points of each support in lexicographicaly decending order

   SortSpt(nVar,nSpt,Spt,Spt1stIdx);

   for (i=0; i<nVar; i++)
   {
      for (j=Spt1stIdx[i]; j<Spt1stIdx[i+1]-1; j++)
      {
         itmp = 0;
         for (k=0; k<nVar; k++) itmp += abs(Spt[j][k]-Spt[j+1][k]);
	 if (itmp == 0)
         {
	    cout << j+1 << "-th and " << j+2 << "-th points in " << i+1 << "-th support are same!" << endl;
	    exit(0);
	 }
      }
   }

   // To locate all non-vertex points in each support
   NonVertex(nVar,nSpt,Spt1stIdx,Spt,NuIdx2OldIdx);

   icnt = -1;
   for (i=0; i<nVar; i++)
   {
      for (j=Spt1stIdx[i]; j<Spt1stIdx[i+1]; j++)
      {
         if (NuIdx2OldIdx[j] == 1)
         {
            icnt = icnt + 1;
            for (k=0; k<nVar; k++) SptVtx[icnt][k] = Spt[j][k];
            NuIdx2OldIdx[icnt] = j;
            //  For the vertex system, NuIdx2OldIdx[i] = j means that 
            //  the current i-th point is the j-th point in the original system
	 }
      }
      SptVtx1stIdx[i+1] = icnt+1;
   }
   SptVtx1stIdx[0] = 0;

   // To reorder the supports if necessary

   int ynReOrder = 0;
   int* OrderSpt = new int [nVar];
   int* iTmp = new int [nVar+1];
   int* iTmp1 = new int [nVar+1];
   int** SptCopy = new int* [Spt1stIdx[nVar]];
   SptCopy[0] = new int [nVar*Spt1stIdx[nVar]];
   for (i=1; i<Spt1stIdx[nVar]; i++) SptCopy[i] = SptCopy[i-1] + nVar;
   int* NuIdx2OldIdxCopy = new int [Spt1stIdx[nVar]];

   for (i=0; i<nVar; i++) iTmp[i] = 0;
   for (i=0; i<nVar; i++)
   {
      itmp = 1000000;
      for (j=0; j<nVar; j++)
      {
         if (iTmp[j] == 0 && SptVtx1stIdx[j+1]-SptVtx1stIdx[j] < itmp)
	 {
            itmp = SptVtx1stIdx[j+1]-SptVtx1stIdx[j];
            OrderSpt[i] = j;
	 }
      }
      if (OrderSpt[i] != i) ynReOrder = 1;
      iTmp[OrderSpt[i]] = 1;
   }

   //  The order of the supports : (OrderSpt[i],i=0,nVar-1)
 
   if (ynReOrder == 1)
   {
      // To change Spt, SptVtx, Spt1stIdx, SptVtx1stIdx
      iTmp[0] = 0;
      iTmp1[0] = 0;
      icnt = -1;
      for (i=0; i<nVar; i++)
      {
         for (j=Spt1stIdx[OrderSpt[i]]; j<Spt1stIdx[OrderSpt[i]+1]; j++)
	 {
	    icnt = icnt + 1;
            for (k=0; k<nVar; k++) SptCopy[icnt][k] = Spt[j][k];
	 }
         iTmp[i+1] = icnt+1;
      }
      for (j=0; j<=icnt; j++)
      {
         for (i=0; i<nVar; i++) Spt[j][i] = SptCopy[j][i];
      }

      icnt = -1;
      for (i=0; i<nVar; i++)
      {
         for (j=SptVtx1stIdx[OrderSpt[i]]; j<SptVtx1stIdx[OrderSpt[i]+1]; j++)
	 {
            icnt = icnt + 1;
            for (k=0; k<nVar; k++) SptCopy[icnt][k] = SptVtx[j][k];
            NuIdx2OldIdxCopy[icnt] = NuIdx2OldIdx[j] - (Spt1stIdx[OrderSpt[i]]-iTmp[i]);
	 }
         iTmp1[i+1] = icnt+1;
      }
      for (j=0; j<=icnt; j++)
      {
         for (i=0; i<nVar; i++) SptVtx[j][i] = SptCopy[j][i];
	 NuIdx2OldIdx[j] = NuIdx2OldIdxCopy[j];
      }

      for (i=0; i<=nVar; i++)  Spt1stIdx[i] = iTmp[i];

      for (i=0; i<=nVar; i++)  SptVtx1stIdx[i] = iTmp1[i];
   }

   // To find whether there are some supports same in the vertex system

   for (i=0; i<nVar; i++) SptType[i] = -1;

   nSpt  = 0;
   for (i=0; i<nVar-1; i++)
   {
      if (SptType[i] < 0)
      {
         nSpt = nSpt + 1;
         for (j=i+1; j<nVar; j++)
	 {
            if (SptType[j] < 0)
	    { 
               if (SptVtx1stIdx[j+1]-SptVtx1stIdx[j] 
	           != SptVtx1stIdx[i+1]-SptVtx1stIdx[i]) continue;
	       int ibreak = 0;
               for (k=SptVtx1stIdx[i]; k<SptVtx1stIdx[i+1]; k++)
	       {
                  for (l=0; l<nVar; l++)
		  {
                     if (SptVtx[k][l] 
		         != SptVtx[k+SptVtx1stIdx[j]-SptVtx1stIdx[i]][l])
		     {
			ibreak = 1;
		        break;
		     }
		  }
	          if (ibreak) break;
	       }
	       if (ibreak) continue;
               SptType[j] = i;
	    }
         }
      }
   }

   if (SptType[nVar-1] < 0) nSpt = nSpt + 1;

   //  StpType[i] = -1 : the i-th support is the base support
   //          = j  : the i-th support is the same as the j-th base support

   if (nSpt<nVar)
   {
      // For the unmixed case, rearrange the supports according to the type
      ReArrangeSpt(nVar,Spt,Spt1stIdx,SptVtx,SptVtx1stIdx,SptType,NuIdx2OldIdx);
   }
   else
   {
      nSpt = nVar;
      for (i=0; i<nVar; i++) SptType[i] = 1;
   }

   delete [] OrderSpt;
   delete [] iTmp;
   delete [] iTmp1;
   delete [] SptCopy[0];
   delete [] SptCopy;
   delete [] NuIdx2OldIdxCopy;
}
#endif
