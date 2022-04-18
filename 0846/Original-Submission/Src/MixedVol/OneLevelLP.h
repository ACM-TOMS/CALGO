#ifndef _OneLevelLP_
#define _OneLevelLP_
 
#include "IT_LP.h"
#include "ExtLP.h"
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>

void OneLevelLP 
   (
   int &Strt1PtTst, int &End1PtTst, bool* PtIn,
   int &CurLPdim, double** A, int &bIdx, double* x, double** Binv, int* Bidx,
   IT_LP &ItLp
   )

{
   int i, j, info;

   extern double* c; //   double* c = new double [CurLPdim];
   extern int* J; //   int* J = new int [n]; // wiht n>=CurLPdim

   // To extend the cells by using 1-point test

   int TstPt = Strt1PtTst;
   for ( i=0; i<CurLPdim; i++ )
      c[i] = -A[TstPt][i];    

   dnulp2_a(End1PtTst+1,CurLPdim,A,bIdx,c,Bidx,x,Binv,info);

   // To record any possilbe points in Bidx passed 1-point test

   j = -1;
   for ( i=0; i<CurLPdim; i++ )
      if ( Bidx[i] >= TstPt )
      {
         PtIn[Bidx[i]] = true;
         J[++j] = Bidx[i];
      }
   if ( ++j > 0 )
   {
      Sort(j,J);
      ItLp.Add(j,J, CurLPdim, Bidx, x, Binv);
   }

   // To do the 1-point test for other points


   if ( info >= 0 )
   {
      for ( TstPt=Strt1PtTst+1; TstPt<=End1PtTst; TstPt++ )
      {
         if ( !PtIn[TstPt] )
         {
            for ( i=0; i<CurLPdim; i++ )
               c[i] = -A[TstPt][i];

            dlp1_1pts(End1PtTst+1,CurLPdim,A,bIdx,c,TstPt,Bidx,x,Binv,
                        PtIn,ItLp);
         }
      }
   }
   else
   {
      for ( TstPt=Strt1PtTst+1; TstPt<=End1PtTst; TstPt++ )
      {
         if ( !PtIn[TstPt] )
         {
            for ( i=0; i<CurLPdim; i++ )
               c[i] = -A[TstPt][i];

            dlp2_1pts(End1PtTst+1,CurLPdim,A,bIdx,c,TstPt,Bidx,x,Binv,
                         PtIn,ItLp);
         }
      }
   }

}

#endif 
