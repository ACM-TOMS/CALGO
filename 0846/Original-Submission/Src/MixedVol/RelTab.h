#ifndef _RelationTable_
#define _RelationTable_  

#include "RlTbLP.h"
#include "L0_IT_LP.h"
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void RelTable
   (
   int &nVar, int &nSpt, int** Spt, int* Spt1stIdx, double* lft,
   bool** RelTab, L0_IML &L0
   )
{
   int    i, j, j1, info, ell;
   double dtmp;
   bool   ynNuPt;
   int    CurSpt, NxSpt, FixPt, TstPt, bIdx, Strt1PtTst, nPtIn;
 
   extern int* J;  //   int* J  = new int [2];
   extern double* c; //   double* c  = new double [nVar+2];

   double** Binv_s  = new double* [nVar];
   Binv_s[0] = new double [ nVar * nVar ];
   for ( i=1; i<nVar; i++ ) Binv_s[i] = Binv_s[i-1] + nVar;
   double* x_s  = new double [nVar];
   int* Bidx_s  = new int [nVar];

   double** a = new double* [Spt1stIdx[nSpt]+1];
   a[0] = new double [ (nVar+2) * (Spt1stIdx[nSpt]+1) ];
   for ( j=1; j<=Spt1stIdx[nSpt]; j++ ) a[j] = a[j-1] + (nVar+2);

   int* LPidx  = new int [Spt1stIdx[nSpt]+1];

   double* x  = new double [nVar+2];
   int* Bidx  = new int [nVar+2];

   double** Binv  = new double* [nVar+2];
   Binv[0] = new double [ (nVar+2) * (nVar+2) ];
   for ( i=1; i<nVar+2; i++ ) Binv[i] = Binv[i-1] + (nVar+2);

   for ( i=0; i<Spt1stIdx[nSpt]; i++ )
      for ( j=0; j<Spt1stIdx[nSpt]; j++ )  RelTab[i][j] = false;

   for ( j=Spt1stIdx[1]; j<Spt1stIdx[nSpt]; j++ ) // [0] changing
   {
      for ( i=0; i<nVar; i++ )  a[j][i] = -Spt[j][i];
      a[j][nVar] = 1.0;   // WRT alpha0
      a[j][nVar+1] = lft[j];
   }
   bIdx = nVar+1;

   for ( CurSpt=0; CurSpt<nSpt; CurSpt++ )
   {
      for ( FixPt=Spt1stIdx[CurSpt]; FixPt<Spt1stIdx[CurSpt+1]; FixPt++ ) 
      {
         nPtIn = -1;            // ! # of constraints involved - 1
         for ( j=Spt1stIdx[CurSpt]; j<FixPt; j++ )
            if ( RelTab[FixPt][j] ) 
            {
               LPidx[++nPtIn] = j;
               for ( i=0; i<nVar; i++ )  a[j][i] = Spt[FixPt][i] - Spt[j][i];
               a[j][nVar] = 1.0;   // for alpha0
               a[j][bIdx] = lft[j] - lft[FixPt];  // the constant term
            }

         Strt1PtTst = nPtIn + 1;    // ! starting index of pts for 1-pt test
         for ( j=FixPt+1; j<Spt1stIdx[CurSpt+1]; j++ ) 
         {
            LPidx[++nPtIn] = j;
            for ( i=0; i<nVar; i++ )  a[j][i] = Spt[FixPt][i] - Spt[j][i];
            a[j][nVar] = 1.0;
            a[j][bIdx] = lft[j] - lft[FixPt];
         }

         // To extend the point by using 1-point test
         // To find a nondegenerate solution for 1-point tests of the points

         RlTbLP2_e(nPtIn+1,nVar,Spt1stIdx[nSpt],a,bIdx,LPidx,Bidx,x,Binv,info );

         if ( FixPt < Spt1stIdx[CurSpt+1]-1 )
         {

            // To record any possilbe points in Bidx passed 1-point test
            if (CurSpt == 0)
            {
               ynNuPt = false;
               for ( i=0; i<nVar; i++ )
                  if ( Bidx[i] > FixPt && !RelTab[FixPt][Bidx[i]] )
                  {
                     ynNuPt = true;
                     RelTab[FixPt][Bidx[i]] = RelTab[Bidx[i]][FixPt] = true;
                  }
               if ( ynNuPt )
               {
                  J[0] = FixPt;
                  for ( i=0; i<nVar; i++ )
                     if ( Bidx[i] > FixPt )
                     {
                        J[1] = Bidx[i];
                        L0.Add(J, nVar, Bidx, x, Binv);   // add 2 points
                     }
               }
            }
            else
               for ( i=0; i<nVar; i++ )
                  if ( Bidx[i] > FixPt ) 
                  {
                     RelTab[Bidx[i]][FixPt] = RelTab[FixPt][Bidx[i]] = true;
                     for ( j=i+1; j<nVar; j++ )
                        if ( Bidx[j] > FixPt )
                           RelTab[Bidx[i]][Bidx[j]] = 
                                 RelTab[Bidx[j]][Bidx[i]] = true;
                  }

            // x, Bidx, Binv are available now
            // To do the 1-point test for other points

            if (info >= 0) 
               for ( j=Strt1PtTst; j<=nPtIn; j++ )
               {
                  TstPt = LPidx[j];
                  if ( !RelTab[FixPt][TstPt] ) 
                  {
                     for ( i=0; i<nVar; i++ )  c[i] = -a[TstPt][i];
                     if (CurSpt == 0) 
                        dlp1_1pt_s(nPtIn+1,nVar,a,bIdx,c,LPidx,FixPt,TstPt,
                                           Bidx,x,Binv,RelTab,L0);
                     else
                        dlp1_1pt_i(nPtIn+1,nVar,a,bIdx,c,LPidx,FixPt,TstPt,
                                           Bidx,x,Binv,RelTab);
                  }
               }
            else
               for ( j=Strt1PtTst; j<=nPtIn; j++ )
               {
                  TstPt = LPidx[j];
                  if ( !RelTab[FixPt][TstPt] )
                  {
                     for ( i=0; i<nVar; i++ )  c[i] = -a[TstPt][i];
                     if ( CurSpt == 0 )
                        dlp2_1pt_s(nPtIn+1,nVar,a,bIdx,c,LPidx,FixPt,TstPt,
                                           Bidx,x,Binv,RelTab,L0);
                     else
                        dlp2_1pt_i(nPtIn+1,nVar,a,bIdx,c,LPidx,FixPt,TstPt,
                                           Bidx,x,Binv,RelTab);
                  }
               }
         }

         for ( i=0; i<nVar; i++ )          // save the starting point for next LP
         {
            for ( j=0; j<nVar; j++ )  Binv_s[i][j] = Binv[i][j];
            Bidx_s[i] = Bidx[i];
            x_s[i] = x[i];
         }

         // To check point FixPt in S_CurSpt and each point in S_i (i>CurSpt)

         nPtIn = -1;  // To get all constraints of pts related to FixPt
         for ( i=Spt1stIdx[CurSpt]; i<Spt1stIdx[CurSpt+1]; i++ )
            if ( RelTab[FixPt][i] )
            {
               LPidx[++nPtIn] = i;
               a[i][nVar] = 0.0;   // no alpha0
            }

         Strt1PtTst = nPtIn + 1;    // starting index of pts for 1-pt test

         for ( NxSpt=CurSpt+1; NxSpt<nSpt; NxSpt++ )
         {
            nPtIn = Strt1PtTst - 1;
            for ( i=Spt1stIdx[NxSpt]; i<Spt1stIdx[NxSpt+1]; i++ )
               LPidx[++nPtIn] = i;
            // The part of Ax<=B for this support was formed at the beginning

            // To do the 1-point test for the points in S_(NxSpt)
            // To form the starting point from saved starting Bidx, x, Binv

            for ( i=0; i<nVar; i++ )
            {
               Bidx[i] = Bidx_s[i];
               x[i] = x_s[i];
            }
            ell = -1;
            x[nVar] = DBL_MAX;
            for ( j1=Strt1PtTst; j1<=nPtIn; j1++ )
            {
               j = LPidx[j1];
               dtmp = a[j][bIdx];
               for ( i=0; i<nVar; i++ )  dtmp -= a[j][i]*x_s[i];
               if ( dtmp < x[nVar] )
               {
                  x[nVar] = dtmp;
                  ell = j;
               }
            }
            Bidx[nVar] = ell;
            for ( i=0; i<nVar; i++ )
            {
               for ( j=0; j<nVar; j++ )  Binv[i][j] = Binv_s[i][j];
               Binv[i][nVar] = 0.0;
            }
            for ( i=0; i<nVar; i++ )  Binv[nVar][i] = 0.0;
            for ( j=0; j<nVar; j++ )
            {
               Binv[j][nVar] = 0.0;
               for ( i=0; i<nVar; i++ )  Binv[j][nVar] -= Binv[j][i]*a[ell][i];
            }
            Binv[nVar][nVar] = 1.0;

            // To find a nondegenerate solution for 1-point tests of the points

            TstPt = LPidx[Strt1PtTst];
            for ( i=0; i<=nVar; i++ )  c[i] = -a[TstPt][i];

            RlTbLP2_a(nPtIn+1,nVar+1,a,bIdx,c,LPidx,Bidx,x,Binv,info);

            // To record any possilbe points in Bidx passed 1-point test

            for ( i=0; i<=nVar; i++ ) 
            {
               if ( Bidx[i] >= Spt1stIdx[NxSpt] )
               {
                  if ( Bidx[i] >= TstPt ) 
                     RelTab[Bidx[i]][FixPt] = RelTab[FixPt][Bidx[i]] = true;
                  for ( j=i+1; j<=nVar; j++ )
                     if ( Bidx[j] >= Spt1stIdx[NxSpt] )
                        RelTab[Bidx[i]][Bidx[j]] 
                                 = RelTab[Bidx[j]][Bidx[i]] = true;
               }
            }

            // To do the 1-point test for other points

            if ( info >= 0) 
               for ( j=Strt1PtTst+1; j<=nPtIn; j++)
               {
                  TstPt = LPidx[j];
                  if ( !RelTab[FixPt][TstPt] ) 
                  {
                     for ( i=0; i<=nVar; i++ )  c[i] = -a[TstPt][i];

                     dlp1_1pt_i(nPtIn+1,nVar+1,a,bIdx,c,LPidx,FixPt,TstPt,
                                        Bidx,x,Binv,RelTab);
                  }
               }
            else
               for ( j=Strt1PtTst+1; j<=nPtIn; j++)
               {
                  TstPt = LPidx[j];
                  if ( !RelTab[FixPt][TstPt] ) 
                  {   
                     for ( i=0; i<=nVar; i++ )  c[i] = -a[TstPt][i];

                     dlp2_1pt_i(nPtIn+1,nVar+1,a,bIdx,c,LPidx,FixPt,TstPt,
                                        Bidx,x,Binv,RelTab);
                  }
               }
         }
      }
   }

   delete [] Binv_s[0];
   delete [] Binv_s;
   delete [] x_s;
   delete [] Bidx_s;

   delete [] a[0];
   delete [] a;
   delete [] x;
   delete [] Bidx;
   delete [] Binv[0];
   delete [] Binv;

   delete [] LPidx;

   return;

}

#endif
