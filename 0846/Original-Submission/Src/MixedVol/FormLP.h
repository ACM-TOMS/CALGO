#ifndef _FORMLP_
#define _FORMLP_
 
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "IT_LP.h"

void FormLP (
            int &nVar,
            int &nSpt, int *SptType, int *Spt1stIdx, bool **RelTab, 
            double** ElmMtrx, int* LstNonZero, int* Lvl2CoDim,
            double ***A, int &bIdx, 
            int &CurLvl, LPdata* ptr, int* CurCell,
            int* Lvl2LPdim, int* Lvl2Spt, 
            bool* Lvl2ynFixFrstPt, bool* Lvl2ynFixLstPt,
            int* MinNumPt, bool** PtIn,
            int &Strt1PtTst, int &End1PtTst, int** LPidx,
            int &CurLPdim, double* x, double** Binv, int* Bidx,
            int &info
            )

{
   bool   ibrk;
   int    ell, i, j, i1, j1, k, kout, EndOldPart, CoDim, FixPt_Orig;
   double smallnum=1.0e-12, dtmp;

   // info : = 0, normal return;
   //        = 1, not enough points to extend, back-tracking.

   // Step 1: To find who will get into the LP constraints

   // CurLvl: index of current level; starting from 1
   // CurSpt: index of current support: starting from 0

   int CurLvl_ = CurLvl - 1;

   CurLPdim = Lvl2LPdim[CurLvl];
   int PreLPdim = Lvl2LPdim[CurLvl-1];

   int CurSpt = Lvl2Spt[CurLvl];

   bool ynFixFrstPt = Lvl2ynFixFrstPt[CurLvl];
   bool ynFixLstPt = Lvl2ynFixLstPt[CurLvl];

   int FixPt = CurCell[CurLvl]; // last point fixed
   if ( CurLvl == 2 )
   {
      FixPt_Orig = FixPt;
      FixPt = Pre2Cur[CurLvl_][FixPt];
   }
   else
      FixPt_Orig = ToOrig[CurLvl_][FixPt];

   int nPtIn = -1;  // # of points involved to form the constraints - 1
   for ( i1=0; i1<FrstPtCurSpt[CurLvl_]; i1++ )
   {
      i = ToOrig[CurLvl_][i1];
      //if ( RelTab[FixPt_Orig][i] && PtIn[CurLvl_][i1] ) 
      if ( RelTab[FixPt_Orig][i] )  // always in
      {
         ToOrig[CurLvl][++nPtIn] = i;
         Cur2Pre[nPtIn] = i1;
         Pre2Cur[CurLvl][i1] = nPtIn;
      }
      else
         Pre2Cur[CurLvl][i1] = -1;   //need??
   }

   FrstPtCurSpt[CurLvl] = nPtIn + 1;
   for ( i1=FrstPtCurSpt[CurLvl_]; i1<FixPt; i1++ )
   {
      i = ToOrig[CurLvl_][i1];
      if ( RelTab[FixPt_Orig][i] && PtIn[CurLvl_][i1] )
      {
         ToOrig[CurLvl][++nPtIn] = i;
         Cur2Pre[nPtIn] = i1;
         Pre2Cur[CurLvl][i1] = nPtIn;
      }
      else
         Pre2Cur[CurLvl][i1] = -1;
   }

   Strt1PtTst = nPtIn + 1;

   for ( i1=FixPt+1; i1<=LstPtCurSpt[CurLvl_]; i1++ )
   {
      i = ToOrig[CurLvl_][i1];
      if ( RelTab[FixPt_Orig][i] && PtIn[CurLvl_][i1] )
      {
         ToOrig[CurLvl][++nPtIn] = i;
         Cur2Pre[nPtIn] = i1;
         Pre2Cur[CurLvl][i1] = nPtIn;
      }
      else
         Pre2Cur[CurLvl][i1] = -1;
   }

   if ( !ynFixLstPt && nPtIn-Strt1PtTst < MinNumPt[CurLvl]-1 )
      { info = 1; return; }  // not enough pts to extend;

   End1PtTst = EndOldPart = LstPtCurSpt[CurLvl] = nPtIn;

   if ( CurLvl == 2 )
      CurCell_Orig[2] = CurCell[2];
   else
      CurCell_Orig[CurLvl] = ToOrig[CurLvl-1][CurCell[CurLvl]];

   if ( ynFixLstPt )    // The last pt fixed is last pt needed;
                        // Extend to next spt.
   {
      Strt1PtTst = FrstPtCurSpt[CurLvl] = nPtIn + 1;

      for ( i=Spt1stIdx[CurSpt+1]; i<Spt1stIdx[CurSpt+2]; i++ )
      {
         ibrk = false;
         for ( j=1; j<=CurLvl; j++ )
         {
            //if ( !RelTab[i][J[j]] )  { ibrk = true; break;  }
            if ( !RelTab[i][CurCell_Orig[j]] )  { ibrk = true; break;  }
         }

         if ( !ibrk )
         {
            ToOrig[CurLvl][++nPtIn] = i;
            Cur2Pre[nPtIn] = i;
            Pre2Cur[CurLvl][i] = nPtIn;
         }
         else
            Pre2Cur[CurLvl][i] = -1;
      }
      End1PtTst = LstPtCurSpt[CurLvl] = nPtIn;

      if ( End1PtTst-Strt1PtTst < SptType[CurSpt+1] )   
         { info = 1; return; }  // not enough pts; need irep(j)+1 pts
   }

   // Step 2: To eliminate a variable by the eq. of the last fixed point

   // To form the equation used for elimination

   double* Elm1Var = v;
   for (i=0; i<PreLPdim; i++ )  Elm1Var[i] = A[CurLvl_][FixPt][i];
   Elm1Var[bIdx] = A[CurLvl_][FixPt][bIdx];

   // To eliminate the kout-th varaible and shift columns of a to left

   if ( !ynFixFrstPt )    // still in same support
   {
      ibrk = false;
      for ( kout=PreLPdim-1; kout>=0; kout-- )
         if ( fabs(Elm1Var[kout]) > smallnum )  { ibrk=true;  break; }

      if ( !ibrk ) 
      {
         cerr << "FormLP.cpp: cannot find a varaible to eliminate.\a\n";
         abort();
      }

      for ( i=0; i<kout; i++ )  Elm1Var[i] /= Elm1Var[kout];
      Elm1Var[bIdx] /= Elm1Var[kout];
//    Elm1Var[kout] = 1.0;

      CoDim = Lvl2CoDim[CurLvl];
      LstNonZero[CoDim] = kout;
      for ( i=0; i<kout; i++ ) ElmMtrx[CoDim][i] = Elm1Var[i];
      ElmMtrx[CoDim][bIdx] = Elm1Var[bIdx];

      for ( i=0; i<=EndOldPart; i++ )
      {
         j = Cur2Pre[i];
         if ( fabs(A[CurLvl_][j][kout]) > smallnum)  // to eliminate
         {
            for ( j1=0; j1<kout; j1++ )
               A[CurLvl][i][j1] = A[CurLvl_][j][j1] 
                                  - A[CurLvl_][j][kout]*Elm1Var[j1];
            A[CurLvl][i][bIdx] = A[CurLvl_][j][bIdx] 
                                  - A[CurLvl_][j][kout]*Elm1Var[bIdx];
         } 
         else
         {
            for ( j1=0; j1<kout; j1++ )  A[CurLvl][i][j1] = A[CurLvl_][j][j1];
            A[CurLvl][i][bIdx] = A[CurLvl_][j][bIdx];
         }
 
         for ( j1=kout; j1<PreLPdim-1; j1++ )  // to shift;
            A[CurLvl][i][j1] = A[CurLvl_][j][j1+1];
      }

      if ( ynFixLstPt )  // set the values for the variable alpha0
      {
         for ( i=0; i<=EndOldPart; i++ )   // constraints without alpha0
            A[CurLvl][i][PreLPdim-1] = 0.0;

         for ( i=FrstPtCurSpt[CurLvl]; i<=nPtIn; i++ ) // constraints with alpha0
         {
            j = ToOrig[CurLvl][i];
            for ( k=0; k<nVar; k++ )
               A[CurLvl][i][k] = A[0][j][k];
            A[CurLvl][i][bIdx] = A[0][j][bIdx];
         }

         for ( i1=1; i1<=CoDim; i1++ )
         {
            kout = LstNonZero[i1];
            Elm1Var = ElmMtrx[i1];

            for ( i=FrstPtCurSpt[CurLvl]; i<=nPtIn; i++ )
            {
               if ( fabs(A[CurLvl][i][kout]) > smallnum)  // to eliminate
               {
                  for ( j1=0; j1<kout; j1++ )
                     A[CurLvl][i][j1] -= A[CurLvl][i][kout]*Elm1Var[j1];
                  A[CurLvl][i][bIdx] -= A[CurLvl][i][kout]*Elm1Var[bIdx];
               }
               for ( j1=kout; j1<nVar-i1; j1++ )  // to shift;
                  A[CurLvl][i][j1] = A[CurLvl][i][j1+1];
            }
         }  // kout = LstNonZero[CoDim] will be used below!

         for ( i=FrstPtCurSpt[CurLvl]; i<=End1PtTst; i++ ) // constraints with alpha0
            A[CurLvl][i][PreLPdim-1] = 1.0;   // For added variable alpha0
      }
   }
   else  // fixing the first point from a support, eliminate alpha0
   {
      kout = PreLPdim - 1;
      if ( fabs(Elm1Var[kout]) > smallnum )  // eliminate alpha0
      {
         for ( i=0; i<kout; i++ )  Elm1Var[i] /= Elm1Var[kout];
         Elm1Var[bIdx] /= Elm1Var[kout];
//       Elm1Var[kout] = 1.0;

         for ( i=0; i<FrstPtCurSpt[CurLvl]; i++ )   // copy the previous a, b saved
         {
            j = Cur2Pre[i];
            for ( j1=0; j1<PreLPdim-1; j1++ )
               A[CurLvl][i][j1] = A[CurLvl_][j][j1];
            A[CurLvl][i][bIdx] = A[CurLvl_][j][bIdx];
         }

         for ( i=FrstPtCurSpt[CurLvl]; i<=End1PtTst; i++ ) // eliminate alpha0 
         {
            j = Cur2Pre[i];
            for ( j1=0; j1<PreLPdim-1; j1++ )
               A[CurLvl][i][j1] = A[CurLvl_][j][j1] - Elm1Var[j1];
            A[CurLvl][i][bIdx] = A[CurLvl_][j][bIdx] - Elm1Var[bIdx];
         }
      }
   }

   ibrk = false;
   if ( CurLvl == 2 )
      for ( i=0; i<PreLPdim; i++ )
      {
         Bidx[i] = ptr->JJJ[i];
         if ( Bidx[i] == FixPt_Orig )
            { k = i;  ibrk = true;  break; }
      }
   else
      for ( i=0; i<PreLPdim; i++ )
      {
         Bidx[i] = ptr->JJJ[i];
         if ( Bidx[i] == FixPt )
            { k = i;  ibrk = true;  break; }
      }

   if ( !ibrk )
   {
      cerr << "FormLP.cpp: no index match for reused info\a\n";
      abort();
   }

   for ( i =k+1; i<PreLPdim; i++ )
   {   
      Bidx[i-1] = ptr->JJJ[i];
   }

   if ( CurLvl > 2 )
   {
      for ( i=0; i<PreLPdim-1; i++ )
         if ( Bidx[i]>-1 ) 
            Bidx[i] = Pre2Cur[CurLvl][Bidx[i]];  // may = -1 if eliminated
   }
   else // = 2
   {
      for ( i=0; i<PreLPdim-1; i++ )
         if ( Bidx[i]>-1 ) 
         {
            Bidx[i] = Pre2Cur[1][Bidx[i]]; // from level 0 to levle 1
            if ( Bidx[i]>-1 ) 
               Bidx[i] = Pre2Cur[CurLvl][Bidx[i]]; // level 1 to level 2
         }
   }

   for ( i=0; i<kout; i++ )  x[i] = ptr->xxx[i];
   for ( i=kout+1; i<PreLPdim; i++ )  x[i-1] = ptr->xxx[i];

   for ( j=0; j<k; j++ )
   {
      for ( i=0; i<kout; i++ )
      {
         Binv[j][i] = ptr->INV[j][i]; 
      }
      for ( i=kout+1; i<PreLPdim; i++ )
      {
         Binv[j][i-1] = ptr->INV[j][i];
      }
   }

   for ( j=k+1; j<PreLPdim; j++ )
   {
      for ( i=0; i<kout; i++ )  Binv[j-1][i] = ptr->INV[j][i]; 
      for ( i=kout+1; i<PreLPdim; i++ )  Binv[j-1][i-1] = ptr->INV[j][i];
   }

   if ( ynFixLstPt )   // the 1st round 1-pt test for next supp
   {
      x[CurLPdim-1] = DBL_MAX;  // The variable alpha0
      ell = -1;
      for ( i=Strt1PtTst; i<=End1PtTst; i++ )
      {
         dtmp = A[CurLvl][i][bIdx];
         for ( i1=0; i1<CurLPdim-1; i1++ )  dtmp -= A[CurLvl][i][i1]*x[i1];

         if (x[CurLPdim-1] > dtmp) 
         {
            ell = i;
            x[CurLPdim-1] = dtmp;
         }
      }
      if ( ell < 0 )
      {
         cerr << "FormLP.cpp: no ell\a\n";
         abort();
      }

      Bidx[CurLPdim-1] = ell;   // constraint becomes the last element of base
      for ( i=0; i<CurLPdim; i++ )               //  update the inverse matrix
         Binv[CurLPdim-1][i] = 0.0;

      for ( j=0; j<CurLPdim-1; j++ )   // -1??? here and next
      {
         dtmp = A[CurLvl][ell][0]*Binv[j][0];
         for ( i=1; i<CurLPdim-1; i++ )  dtmp += A[CurLvl][ell][i]*Binv[j][i];

         Binv[j][CurLPdim-1] = -dtmp;
      }
      Binv[CurLPdim-1][CurLPdim-1] = 1.0;
   }

   for ( i=FrstPtCurSpt[CurLvl]; i<Strt1PtTst; i++ )
      PtIn[CurLvl][i] = true;
   for ( i=Strt1PtTst; i<=End1PtTst; i++ )
      PtIn[CurLvl][i] = false;

   info = 0;
   return;
}

#endif
