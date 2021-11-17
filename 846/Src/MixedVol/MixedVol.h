#ifndef _Mixed_Vol_
#define _Mixed_Vol_

// Global variables

int CurLvl;

int* J;   // for ItLp.Add

double* v; // for LPs
double* c;

int** iwork;  // for CellVol

int* CurCell_Orig;
int* FrstPtCurSpt;
int* LstPtCurSpt;
int** ToOrig;
int* Cur2Pre;
int** Pre2Cur;

// End of global variables

#include "Sort.h"
#include "FormLP.h"
#include "OneLevelLP.h"
#include "CellVol.h"
#include "CellStack.h"
#include "L0_IT_LP.h"
#include "IT_LP.h"
#include "RelTab.h"
#include "FormLP1.h"


void MixedVol 
   (
   int &nVar, int &nSpt, int* SptType, int* Spt1stIdx,
   int** Spt, double* lft, CellStack &MCells, int &MVol
   )
              
{
   MVol=0;

   int i, j, k; 
   int CurLPdim, info, Strt1PtTst, End1PtTst;

   v = new double [nVar+1];
   c = new double [nVar+2];
   iwork = new int* [nVar];
   iwork[0] = new int [nVar*nVar];
   for (int i=1; i<nVar; i++) iwork[i] = iwork[i-1] + nVar;
   int* iwork_head = iwork[0];

   int nMCells=0;

   int CellSize = nSpt;
   for ( i=0; i<nSpt; i++ )
      CellSize += SptType[i];

   IT_LP      ItLp( nSpt, SptType );
   L0_IML     L0;

   J  = new int [CellSize];

   bool** RelTab  = new bool* [Spt1stIdx[nSpt]];
   RelTab[0] = new bool [ Spt1stIdx[nSpt] * Spt1stIdx[nSpt] ];
   for ( i=1; i<Spt1stIdx[nSpt]; i++ ) 
      RelTab[i] = RelTab[i-1] + Spt1stIdx[nSpt];

   RelTable(nVar, nSpt, Spt, Spt1stIdx, lft, RelTab, L0);

   LPdata* ptr;

   k = 0;
   J[0] = Spt1stIdx[nSpt];
   info = Spt1stIdx[nSpt]; // # of rows in A
   for ( i=0; i<nSpt-1; i++ )
   {
      for ( j=0; j<SptType[i]; j++ )
      {
         J[++k] = Spt1stIdx[i+1];
         info += Spt1stIdx[i+1];
      }
      J[++k] = Spt1stIdx[i+2];
      info += Spt1stIdx[i+2];
   }
   for ( j=0; j<SptType[nSpt-1]; j++ )
   {
      J[++k] = Spt1stIdx[nSpt];
      info += Spt1stIdx[nSpt];
   }

   double*** A = new double** [CellSize];
   A[0] = new double* [ info ];
   A[0][0] = new double [ info * (nVar+1) ];
   info = 0;
   for ( j=1; j<J[0]; j++ )
      A[0][j] = A[0][j-1] + (nVar+1);
   for ( i=1; i<CellSize; i++ )
   {
      A[i] = A[i-1] + J[i-1];
      A[i][0] = A[i-1][J[i-1]-1] + (nVar+1);
      for ( j=1; j<J[i]; j++ )
         A[i][j] = A[i][j-1] + (nVar+1);
   }

   int** LPidx  = new int* [CellSize];
   LPidx[0]  = new int [ CellSize * Spt1stIdx[nSpt] ];
   for ( i=1; i<CellSize; i++ ) LPidx[i] = LPidx[i-1] + Spt1stIdx[nSpt];
   bool** PtIn  = new bool* [CellSize];
   PtIn[0]  = new bool [ CellSize * Spt1stIdx[nSpt] ];
   for ( i=1; i<CellSize; i++ ) PtIn[i] = PtIn[i-1] + Spt1stIdx[nSpt];

   double* x  = new double [nVar];
   int* Bidx  = new int [nVar];
   double** Binv  = new double* [nVar];
   Binv[0] = new double [ nVar * nVar ];
   for ( i=1; i<nVar; i++ ) Binv[i] = Binv[i-1] + nVar;

   double** ElmMtrx  = new double* [CellSize];
   ElmMtrx[0] = new double [ CellSize * (nVar+1) ];
   for ( i=1; i<CellSize; i++ ) ElmMtrx[i] = ElmMtrx[i-1] + (nVar+1);
   int* LstNonZero = new int [CellSize];
   int* Lvl2CoDim  = new int [CellSize];
   int* CurCell  = new int [CellSize];
   CurCell_Orig  = new int [CellSize];
   int* Lvl2LPdim  = new int [CellSize];
   int* Lvl2Spt  = new int [CellSize];
   int* Lvl2MinNumPt  = new int [CellSize];
   bool* Lvl2ynFixFrstPt  = new bool [CellSize];
   bool* Lvl2ynFixLstPt  = new bool [CellSize];
   FrstPtCurSpt = new int [CellSize];
   LstPtCurSpt = new int [CellSize];
   ToOrig  = new int* [CellSize];
   ToOrig[0]  = new int [ CellSize * Spt1stIdx[nSpt] ];
   for ( i=1; i<CellSize; i++ ) ToOrig[i] = ToOrig[i-1] + Spt1stIdx[nSpt];
   Cur2Pre  = new int [ Spt1stIdx[nSpt] ];
   Pre2Cur  = new int* [CellSize];
   Pre2Cur[0]  = new int [ CellSize * Spt1stIdx[nSpt] ];
   for ( i=1; i<CellSize; i++ ) Pre2Cur[i] = Pre2Cur[i-1] + Spt1stIdx[nSpt];

   // Note that the level 0 is the artificial head
   Lvl2LPdim[0] = nVar + 1; //just for conveniece, should = nVar, see below
   Lvl2CoDim[0] = 0; // # of variables eliminated. [0] not used
   k = 0;
   for ( i=0; i<nSpt; i++ )
   {
      ++k;
      Lvl2LPdim[k] = Lvl2LPdim[k-1] - 1;
      Lvl2CoDim[k] = Lvl2CoDim[k-1];
      Lvl2Spt[k] = i;
      Lvl2MinNumPt[k] = SptType[i];
      Lvl2ynFixFrstPt[k] = true;
      Lvl2ynFixLstPt[k] = false;
      for ( j=1; j<=SptType[i]; j++ )
      {
         ++k;
	 if ( k==CellSize ) break;
         Lvl2LPdim[k] = Lvl2LPdim[k-1] - 1;
         Lvl2CoDim[k] = Lvl2CoDim[k-1] + 1;
         Lvl2Spt[k] = i;
         Lvl2MinNumPt[k] = Lvl2MinNumPt[k-1] -1;
         Lvl2ynFixFrstPt[k] = false;
         Lvl2ynFixLstPt[k] = false;
      }
      if ( k==CellSize ) break;
      Lvl2LPdim[k] += 1;  // add alpha0 for next spt
      if ( i < nSpt-1 ) Lvl2MinNumPt[k] = SptType[i+1] + 1;
      Lvl2ynFixLstPt[k] = true;
   }      
   Lvl2LPdim[0] = nVar; //=nVar+1-1, in RelTab, add alpha0, fix one point

   //To define A[0]
   for ( i=Spt1stIdx[0]; i<Spt1stIdx[nSpt]; i++ )
   {
      for ( j=0; j<nVar; j++ ) A[0][i][j] = -Spt[i][j];
      A[0][i][nVar] = lft[i];  //bIdx = nVar
   }
   //To define PtIn[0]
   for ( i=Spt1stIdx[0]; i<Spt1stIdx[1]; i++ )
      PtIn[0][i] = true;

   while(  L0.Migrate(ItLp.ResetCurLevelTo1()) )
   {
      ItLp.RenewNP1();

      while( ! ItLp.IsEmpty() )
      {
         if( ! ItLp.IsFull() )
         {
            CurLvl = ItLp.Level();  // note that level-0 is the artificial head
            CurCell[CurLvl] = ItLp.FixedIdxNdPtr()->idx; 
                                    // The index of the point to be fixed

            ptr = ItLp.FixedIdxNdPtr()->info; // The ptr to saved x, Binv and jj

            if ( CurLvl <= 1 )
            {
               // To do elimination and form new system Ax<=b

               CurCell_Orig[CurLvl] = CurCell[CurLvl];

               FormLP1( nSpt,SptType,Spt1stIdx,RelTab,
                        A,nVar,
                        CurLvl,CurCell,Lvl2LPdim,
                        Lvl2Spt,Lvl2ynFixLstPt,
                        Lvl2MinNumPt, PtIn,
                        LPidx,
                        info);
            }
            else
            {
               // To do elimination, form new system Ax<=b and then solve the LP

               FormLP(  nVar,
                        nSpt,SptType,Spt1stIdx,RelTab,
                        ElmMtrx,LstNonZero,Lvl2CoDim,
                        A,nVar,
                        CurLvl,ptr,CurCell,Lvl2LPdim,Lvl2Spt,
                        Lvl2ynFixFrstPt, Lvl2ynFixLstPt,
                        Lvl2MinNumPt, PtIn,
                        Strt1PtTst,End1PtTst,LPidx,
                        CurLPdim, x, Binv, Bidx,
                        info);

               if (info == 0 )
                  OneLevelLP(Strt1PtTst,End1PtTst,PtIn[CurLvl],
                              CurLPdim, A[CurLvl], nVar, x, Binv, Bidx,
                              ItLp);
            }
         } // This level finished
         while( !ItLp.NextLevel() && !ItLp.IsEmpty() )
         {
            if( ItLp.IsFull() )
            {
               // To store the cell and calculate its volume
               J[0] = ItLp.Cell()[0];
               J[1] = ItLp.Cell()[1];
               for ( i=2; i<CellSize; i++ ) J[i] = ToOrig[i][ItLp.Cell()[i]];
               MCells.Push( J ); //Stores the cell (indices of the vertex support)

               CellVol(nVar,nSpt,Spt,SptType,J,MVol);
               ++nMCells;  // nMCells = # of mixed cells at the end
            }
            ItLp.StepBack();
         }
      } // End of while ( !ItLp.IsEmpty() )
   } // End of while ( L0.Migrate(inp) )

   // clear memory space
   delete [] Binv[0];
   delete [] Binv;
   delete [] x;
   delete [] RelTab[0];
   delete [] RelTab;

   delete [] A[0][0];
   delete [] A[0];
   delete [] A;

   delete [] LPidx[0];
   delete [] LPidx;
   delete [] PtIn[0];
   delete [] PtIn;
   delete [] Bidx;
   delete [] J;
   delete [] v;
   delete [] c;
 
   delete [] CurCell;
   delete [] CurCell_Orig;
   delete [] Lvl2LPdim;
   delete [] Lvl2Spt;
   delete [] Lvl2MinNumPt;
   delete [] Lvl2ynFixFrstPt;
   delete [] Lvl2ynFixLstPt;
   delete [] ElmMtrx[0];
   delete [] ElmMtrx;
   delete [] LstNonZero;
   delete [] Lvl2CoDim;

   delete [] FrstPtCurSpt;
   delete [] LstPtCurSpt;
   delete [] ToOrig[0];
   delete [] ToOrig;
   delete [] Cur2Pre;
   delete [] Pre2Cur[0];
   delete [] Pre2Cur;

   delete [] iwork_head;
   delete [] iwork;

   return;
}
#endif
