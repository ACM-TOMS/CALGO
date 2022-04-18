#ifndef _FORMLP1_
#define _FORMLP1_
 
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>

void FormLP1(
            int &nSpt, int *SptType, int *Spt1stIdx, bool **RelTab, 
            double ***A, int &bIdx, 
            int &CurLvl, int* CurCell,
            int* Lvl2LPdim, int* Lvl2Spt,
            bool* Lvl2ynFixLstPt,
            int* MinNumPt, bool** PtIn,
            int** LPidx,
            int &info
            )

{
   int    i, j, j1;

   // info : = 0, normal return;
   //        = 1, not enough points to extend, back-tracking.

   // Step 1: To find who will get into the LP constraints

   // CurLvl: index of current level; starting from 1 ???

   int CurLvl_ = CurLvl - 1;

   int PreLPdim = Lvl2LPdim[CurLvl-1];

   bool ynFixLstPt = Lvl2ynFixLstPt[CurLvl];

   int FixPt = CurCell[CurLvl]; // last point fixed

   int nPtIn = -1;  // # of points involved to form the constraints - 1
   FrstPtCurSpt[CurLvl] = 0;

   for ( i=Spt1stIdx[0]; i<FixPt; i++ )
      if ( RelTab[FixPt][i] )
      {
         ToOrig[CurLvl][++nPtIn] = i;
         Pre2Cur[CurLvl][i] = nPtIn;  // need for saved Bidx for next lvl
      }
      else
         Pre2Cur[CurLvl][i] = -1;  // need??
   Pre2Cur[CurLvl][FixPt] = -1;

   int Strt1PtTst = nPtIn + 1;

   for ( i=FixPt+1; i<Spt1stIdx[1]; i++ )
      if ( RelTab[FixPt][i] )
      {
         ToOrig[CurLvl][++nPtIn] = i;
         Pre2Cur[CurLvl][i] = nPtIn;   // need for saved Bidx for next lvl
      }
      else
         Pre2Cur[CurLvl][i] = -1;

   if ( !ynFixLstPt && nPtIn-Strt1PtTst < MinNumPt[CurLvl]-1 )
      { info = 1; return; }  // not enough pts to extend;

   LstPtCurSpt[CurLvl] = nPtIn;

   // To form the matrix

   for ( i=0; i<=LstPtCurSpt[CurLvl]; i++ )   // copy the previous a, b saved
   {
      j = ToOrig[CurLvl][i];
      for ( j1=0; j1<PreLPdim; j1++ )  
         A[CurLvl][i][j1] = A[CurLvl_][j][j1] - A[CurLvl_][FixPt][j1];

      A[CurLvl][i][bIdx] = A[CurLvl_][j][bIdx] - A[CurLvl_][FixPt][bIdx];
   }

   for ( i=0; i<=LstPtCurSpt[CurLvl]; i++ )
      PtIn[CurLvl][i] = true;

   info = 0;
   return;
}

#endif
