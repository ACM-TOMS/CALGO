#ifndef _Relation_Table_LPCODE_
#define _Relation_Table_LPCODE_  

#include <string.h>
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "IT_LP.h"
#include "L0_IT_LP.h"

void RlTbLP2_a 
   (
   int ma, int na, double** a, int &bIdx, double* c,
   int* LPidx, int* Bidx, double* x, double** Binv, 
   int &info
   )
{

   bool   ibrk;
   int    ell, k, k_, i, i1, j, flag;
   double eps=1.0e-6, sigj, vmax, bot, top, dtmp;

   extern double* v; //   double* v = new double [na];   
 
   while(1)
   { 
      // step 1: To find which constraint to kick out
      ibrk = false;
      for ( i=0; i<na; i++ ) if ( Bidx[i] == -1 ) { ibrk = true;  break; }

      if ( ibrk )
      {
         vmax = -1.0;
         for ( i=0; i<na; i++ )
            if ( Bidx[i] == -1 )
            {
               dtmp = c[0]*Binv[i][0];
               for ( j=1; j<na; j++ ) dtmp += c[j]*Binv[i][j];
               if ( fabs(dtmp) > vmax )  { bot = dtmp;  vmax = fabs(bot);  k = i; }
            }
         if ( vmax > eps )
            if ( bot >= 0.0 )
               for ( i=0; i<na; i++ ) v[i] = Binv[k][i];
            else
               for ( i=0; i<na; i++ ) v[i] = -Binv[k][i];
         else
            ibrk = false;
      }
      if ( !ibrk )
      {
         vmax = -DBL_MAX;
         for ( i=0; i<na; i++ )
            if ( Bidx[i] != -1 )
            {
               dtmp = c[0]*Binv[i][0];
               for ( j=1; j<na; j++ ) dtmp += c[j]*Binv[i][j];
               if ( dtmp > vmax ) { vmax = dtmp;  k = i; }
            }
         if ( vmax < eps ) break;  // found an optimal solution
                                    // go to solve the smaller LP

         for ( i=0; i<na; i++ )  v[i] = Binv[k][i];
      }

      // step 2: To find which constraint coming in

      sigj = DBL_MAX;
      ell = -1;
      for ( i1=0; i1<ma; i1++ )
      {
         i = LPidx[i1];
         ibrk = false;
         for ( j=0; j<na; j++ )  if ( Bidx[j] == i )  { ibrk = true;  break; }
         if ( ibrk )  continue;

         bot = a[i][0]*v[0];
         for ( j=1; j<na; j++ ) bot += a[i][j]*v[j];
         if ( bot >= -eps )  continue;

         top = -a[i][bIdx];
         for ( j=0; j<na; j++ )  top += a[i][j]*x[j];
         top /= bot;
         if ( top >= sigj )  continue;

         sigj = top;
         ell = i;
      }
 
      if ( ell < 0 )
      {
         cerr << "dnulp2_a.cpp: LP unbounded\n\a";
         abort();
      }

      // step 3: To update x, Bidx, Binv

      for ( i=0; i<na; i++ )  x[i] -= sigj * v[i];
      top = a[ell][0] * Binv[k][0];
      for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[k][i];
      if ( fabs(top) < eps )
      {
         cerr << "dnulp2_a.cpp: Base matrix is singular\a\n";
         abort();
      }

      top = 1.0 / top;
      for ( i=0; i<na; i++ )  Binv[k][i] *= top;
      for ( j=0; j<na; j++ )
         if ( j != k )
         {
            top = a[ell][0] * Binv[j][0];
            for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[j][i] ;
            for ( i=0; i<na; i++ )  Binv[j][i] -= top * Binv[k][i];
         }

      Bidx[k] = ell;
   }

   k_ = -1;
   info = 0;
   for ( i=0; i<na; i++ )
      if ( Bidx[i] > -1 )  ++info;

   if ( info == na ) return;   // case of full rank

   k = 0;
   while ( k < na )
   {
      if ( Bidx[k] > -1 || k == k_ )  { ++k;  continue; }

      for ( i=0; i<na; i++ ) v[i] = Binv[k][i];

      flag = 1;

      // step 2: To find which constraint coming in

      ell = -1;
      while ( ell == -1 )
      {
         sigj = DBL_MAX;
         for ( i1=0; i1<ma; i1++ )
         {
            i = LPidx[i1];
            ibrk = false;
            for ( j=0; j<na; j++ )  if ( Bidx[j] == i )  { ibrk = true;  break; }
            if ( ibrk )  continue;

            bot = a[i][0]*v[0];
            for ( j=1; j<na; j++ )  bot += a[i][j]*v[j];
            if ( bot >= -eps )  continue;

            top = -a[i][bIdx];
            for ( j=0; j<na; j++ )  top += a[i][j]*x[j];
            top /= bot;
            if ( top >= sigj )  continue;

            sigj = top;
            ell = i;
         }

         ibrk = false;
         if ( ell < 0 )
            if ( flag == 1 )
            {
               for ( i=0; i<na; i++ ) v[i] = -Binv[k][i];
               flag = 0;
            }
            else
            {
               k = k + 1;  ibrk = true;  break;
            }
      }
 
      if ( ibrk ) continue;

      // step 3:  To update x, Bidx, Binv

      for ( i=0; i<na; i++ )  x[i] -= sigj * v[i];

      top = a[ell][0] * Binv[k][0];
      for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[k][i];
      if( fabs(top) < eps )
      {
         cerr << "dnulp2_a.cpp: Base matrix is singular\a\n";
         abort();
      }
      top = 1.0 / top;
      for ( i=0; i<na; i++ )  Binv[k][i] *= top;
      for ( j=0; j<na; j++ )
         if ( j != k )
         {
            top = a[ell][0] * Binv[j][0];
            for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[j][i] ;
            for ( i=0; i<na; i++ )  Binv[j][i] -= top * Binv[k][i];
         }
         Bidx[k] = ell;
         info = info + 1;

         if ( info == na ) 
            return;
         else
         { k_ = k;  k = 0; }   // necessary?
   }
   info = - info;   // the case of not full rank
}

////////////////////////////////////////

#include <string.h>
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "IT_LP.h"
#include "L0_IT_LP.h"

void RlTbLP2_e
   (
   int ma_, int na_, int &NumCol, double** a, 
   int &bIdx, int* LPidx, int* Bidx, double* x, double** Binv, 
   int &info
   )
{
   bool   ibrk;
   int    ell, k, k_, i, i1, j, na, ma, flag;
   double eps=1.0e-6, sigj, vmax, bot, top, dtmp;

   extern double* v; //   double* v = new double [na];

   // Expand the system (add variable epsilon)
   // Assumed a, b, Binv, x have extra 1 dimension to hold epsilon

   na = na_ + 1;
   ma = ma_ + 1;
   for ( i=0; i<na_; i++ )  a[NumCol][i] = 0.0;  // -epsilon<=0
   a[NumCol][na_] = -1.0;
   a[NumCol][bIdx] = 0.0;
   LPidx[ma_] = NumCol;
   dtmp = 0.0;
   ell = NumCol;             // the index of the equality equation
   for ( i=0; i<ma_; i++ )
   {
      j = LPidx[i];
      if ( a[j][bIdx] < -eps )
      {
         a[j][na-1] = -1;
         if ( -a[j][bIdx] > dtmp ) { dtmp = -a[j][bIdx];  ell = j; }
      }
      else
         a[j][na-1] = 0.0;
   } 

   // To form the first x, Bidx, Binv
   for ( i=0; i<na_; i++ ) Bidx[i] = -1;  Bidx[na-1] = ell;
   for ( i=0; i<na_; i++ ) x[i] = 0.0;  x[na_] = dtmp;
   for (i=0; i<na; i++)
   {
      for (j=0; j<na; j++)  Binv[i][j] = 0.0;
      Binv[i][i] = 1.0;
   }
   for (j=0; j<na_; j++) Binv[j][na-1] = a[ell][j];
   Binv[na-1][na-1] = -1.0;

   // To apply the LP algorithm to the larger LP
   // step 1: To find which constraint to kick out
 

   while(1)
   {
      ibrk = false;
      for ( i=0; i<na; i++ ) if ( Bidx[i] == -1 ) { ibrk = true;  break; }

      if ( ibrk )
      {
         vmax = -1.0;
         for ( i=0; i<na; i++ )
            if ( Bidx[i] == -1 ) 
            {
               dtmp = Binv[i][na-1];
               if ( fabs(dtmp) > vmax )
               {
                  bot = dtmp;  vmax = fabs(bot);  k = i;
               }
            }
         if ( vmax > eps )
            if ( bot >= 0.0 )
               for ( i=0; i<na; i++ ) v[i] = Binv[k][i];
            else
               for ( i=0; i<na; i++ ) v[i] = -Binv[k][i];
         else
            ibrk = false;
      }
      if ( !ibrk )
      {
         vmax = -DBL_MAX;
         for ( i=0; i<na; i++ )
            if ( Bidx[i] != -1 ) 
            {
               dtmp = Binv[i][na-1];
               if (dtmp > vmax) { vmax = dtmp;  k = i; }
            }
         if (vmax < eps) break;  // found an optimal solution
                                 // go to solve the smaller LP
         for ( i=0; i<na; i++ ) v[i] = Binv[k][i];
      }

      // step 2: To find which constraint coming in
      sigj = DBL_MAX;
      ell = -1;
      for ( i1=0; i1<ma; i1++)
      {
         i = LPidx[i1];
         ibrk = false;
         for ( j=0; j<na; j++ )  if ( Bidx[j] == i )  {ibrk = true;  break;}
         if ( ibrk )  continue;

         bot = a[i][0]*v[0];
         for ( j=1; j<na; j++ ) bot += a[i][j]*v[j];
         if ( bot >= -eps )  continue;

         top = -a[i][bIdx];
         for ( j=0; j<na; j++ )  top += a[i][j]*x[j];
         top /= bot;
         if ( top >= sigj )  continue;

         sigj = top;
         ell = i;
      }

      if ( ell < 0 )
      {
         cerr << "dnulp2_e.cpp: LP unbounded??\n\a";
         abort();
      }

      // step 3: To update x, Bidx, Binv
      for ( i=0; i<na; i++ )  x[i] -= sigj * v[i];
      top = a[ell][0] * Binv[k][0];
      for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[k][i];
      if ( fabs(top) < eps )
      {
         cerr << "dnulp2_e.cpp: Base matrix is singular\a\n";
         abort();
      }
      top = 1.0 / top;
      for ( i=0; i<na; i++ )  Binv[k][i] *= top;
      for ( j=0; j<na; j++ )
         if ( j != k )
         {
            top = a[ell][0] * Binv[j][0];
            for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[j][i] ;
            for ( i=0; i<na; i++ )  Binv[j][i] -= top * Binv[k][i];
         }
      Bidx[k] = ell;
   }

   // To change x, Bidx, Binv into the ones without the x(na)=epsilon variable
   if ( Bidx[na-1] != NumCol )
   {
      k = -1;
      j = -1;
      while ( j < na-1 )
      {
         j = j + 1;
         if ( Bidx[j] == NumCol )  { k = j;  j = na-1;  break; }
      }
      if ( k == -1 )   // for testing
      {
         cerr << "dnulp2_e.cpp: no index ma in Bidx \a\n";
         abort();       // not happen in our case
      }

      Bidx[k] = Bidx[na-1];
      for ( i=0; i<na_; i++ ) Binv[k][i] = Binv[na-1][i];
   }

   for ( i=0; i<ma_; i++ )  a[LPidx[i]][na-1] = 0.0;  //reset values WRT epsilon

   info = 0;
   for ( i=0; i<na_; i++ )
      if ( Bidx[i] > -1)  ++info;

   if ( info == na_ ) return;   // case of full rank

   k_ = -1;
   k = 0;
   while ( k < na_ )
   { 
      if ( Bidx[k] > -1 || k == k_ )  {  ++k;  continue; }

      for ( i=0; i<na_; i++ )  v[i] = Binv[k][i];

      flag = 1;
 
      // step 2: To find which constraint coming in
      ell = -1;
      while ( ell == -1 )
      {
         sigj = DBL_MAX;
         for ( i1=0; i1<ma_; i1++ )
         {
            i = LPidx[i1];
            ibrk = false;
            for ( j=0; j<na_; j++ ) if ( Bidx[j] == i )  { ibrk = true;  break; }
            if ( ibrk )  continue;

            bot = a[i][0]*v[0];
            for ( j=1; j<na_; j++ ) bot += a[i][j]*v[j];
            if ( bot >= -eps )  continue;

            top = -a[i][bIdx];
            for ( j=0; j<na_; j++ )  top += a[i][j]*x[j];
            top /= bot;
            if ( top >= sigj )  continue;

            sigj = top;
            ell = i;
         }
         ibrk = false;
         if ( ell < 0 ) 
            if ( flag == 1)
            {
               for ( i=0; i<na_; i++ ) v[i] = -Binv[k][i];
               flag = 0;
            }
            else
            {
               k = k + 1;  ibrk = true;  break;
            }
      }

      if ( ibrk ) continue;

      // step 3:  To update x, Bidx, Binv
      for ( i=0; i<na_; i++ )  x[i] -= sigj * v[i];

      top = a[ell][0] * Binv[k][0];
      for ( i=1; i<na_; i++ )  top += a[ell][i] * Binv[k][i];
      if( fabs(top) < eps )
      {
         cerr << "dnulp2_e.cpp: Base matrix is singular\a\n";
         abort();
      } 
      top = 1.0 / top;
      for ( i=0; i<na_; i++ )  Binv[k][i] *= top;
      for ( j=0; j<na_; j++ )
         if ( j != k )
         {
            top = a[ell][0] * Binv[j][0];
            for ( i=1; i <na_; i++ )  top += a[ell][i] * Binv[j][i] ;
            for ( i=0; i <na_; i++ )  Binv[j][i] -= top * Binv[k][i];
         }
      Bidx[k] = ell;
      ++info;

      if ( info == na_ )  
         return;
      else 
         { k_ = k;  k = 0; }   // necessary?
   }
   info = - info;   // the case of not full rank
}

////////////////////////////////////////

#include <string.h>
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "IT_LP.h"
#include "L0_IT_LP.h"

void dlp2_1pt_i
   (
   int ma, int na, double** a, int &bIdx, 
   double* c, int* LPidx, int &FixPt, int &TstPt, 
   int* Bidx, double* x, double** Binv, bool** RelTab
   )

{
   bool ibrk;
   int    i, i1, j, k, ell;
   double eps=1.0e-6, vmax, top, bot, sigj, dtmp ;
 
   while(1) // LP loop
   {
      // STEP 1: FIND WHICH CONSTRAINT BEING KICKED OUT
      vmax = -DBL_MAX;
      for ( i=0; i<na; i++ )
         if ( Bidx[i] != -1 )
         {
            dtmp = c[0]*Binv[i][0];
            for ( j=1; j<na; j++ )  dtmp += c[j]*Binv[i][j];
            if ( dtmp > vmax )  { vmax = dtmp;  k = i; }
         }
 
      if ( vmax < eps ) return;      // out with optimal solution

      // STEP 2: FIND WHICH CONSTRAINT COMING IN
      sigj = DBL_MAX;
      ell = -1;
      for ( i1=0; i1<ma; i1++ )
      {
         i = LPidx[i1];
         ibrk = false;
         for ( j=0; j<na; j++ )  if ( Bidx[j] == i )  { ibrk = true;  break; }
         if ( ibrk )  continue;

         bot = a[i][0]*Binv[k][0];
         for ( j=1; j<na; j++ ) bot += a[i][j]*Binv[k][j];
         if ( bot >= -eps )  continue;

         top = -a[i][bIdx];
         for ( j=0; j<na; j++ )  top += a[i][j]*x[j];
         top /= bot;
         if ( top >= sigj )  continue;

         sigj = top;
         ell = i;
      }

      if ( ell < 0 )
      {
         cerr << "dlp2_1pt_i.cpp: LP unbounded\n\a";
         abort();       
      }

      // STEP 3: UPDATE x, Bidx, Binv
      for ( i=0; i<na; i++ )  x[i] -= sigj * Binv[k][i];

      top = a[ell][0] * Binv[k][0];
      for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[k][i];
      if ( fabs(top) < eps )
      {
         cerr << "dlp2_1pt_i.cpp: Base matrix is singular\a\n";
         abort();
      }
      top = 1.0 / top;
      for ( i=0; i<na; i++ )  Binv[k][i] *= top;
      for ( j=0; j<na; j++ )
         if ( j != k )
         {
            top = a[ell][0] * Binv[j][0];
            for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[j][i];
            for ( i=0; i<na; i++ )  Binv[j][i] -= top * Binv[k][i];
         }

      Bidx[k] = ell;

      if( ell >= TstPt ) 
      {
         for ( i=0; i< na; i++ )
            if(i != k && Bidx[i] > -1 && !RelTab[ell][Bidx[i]] )
               RelTab[ell][Bidx[i]] = RelTab[Bidx[i]][ell] = true;

         if ( !RelTab[ell][FixPt] ) 
            RelTab[ell][FixPt] = RelTab[FixPt][ell] = true;
      }
   } // LP loop
} // EOF


////////////////////////////////////////

#include <string.h>
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "IT_LP.h"
#include "L0_IT_LP.h"

void dlp2_1pt_s
   (
   int ma, int na, double** a, int &bIdx, 
   double* c, int* LPidx, int &FixPt, int &TstPt, 
   int* Bidx, double* x, double** Binv, bool** RelTab,
   L0_IML &L0
   )

{
   bool   ibrk;
   int    i, i1, j, k, ell;
   double eps=1.0e-6, vmax, top, bot, sigj, dtmp;
 
   extern int* J;  //   int* J = new int [2];

   while(1) // LP loop
   {
      // STEP 1: FIND WHICH CONSTRAINT BEING KICKED OUT
      vmax = -DBL_MAX;
      for ( i=0; i<na; i++ ) 
         if ( Bidx[i] != -1 )
         {
            dtmp = c[0]*Binv[i][0];
            for ( j=1; j<na; j++ )  dtmp += c[j]*Binv[i][j];
            if ( dtmp > vmax )  { vmax = dtmp;  k = i; }
         }
 
      if ( vmax < eps ) return;      // out with optimal solution

      // STEP 2: FIND WHICH CONSTRAINT COMING IN
      sigj = DBL_MAX;
      ell = -1;
      for ( i1=0; i1<ma; i1++ )
      {
         i = LPidx[i1];
         ibrk = false;
         for ( j=0; j<na; j++ )  
            if ( Bidx[j] == i )  { ibrk = true;  break; }
         if ( ibrk )  continue;

         bot = a[i][0]*Binv[k][0];
         for ( j=1; j<na; j++ )  bot += a[i][j]*Binv[k][j];
         if ( bot >= -eps )  continue;

         top = -a[i][bIdx];
         for ( j=0; j<na; j++ )  top += a[i][j]*x[j];
         top /= bot;
         if ( top >= sigj )  continue;

         sigj = top;
         ell = i;
      }

      if ( ell < 0 )
      {
         cerr << "dlp2_1pt_s.cpp: LP unbounded\n\a";
         abort();       
      }

      // STEP 3: UPDATE x, Bidx, Binv
      for ( i=0; i<na; i++ )  x[i] -= sigj * Binv[k][i];

      top = a[ell][0] * Binv[k][0];
      for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[k][i];
      if ( fabs(top) < eps )
      {
         cerr << "dlp2_1pt_s.cpp: Base matrix is singular\a\n";
         abort();
      }
      top = 1.0 / top;
      for ( i=0; i<na; i++ )  Binv[k][i] *= top;
      for ( j=0; j<na; j++ )
         if ( j != k )
         {
            top = a[ell][0] * Binv[j][0];
            for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[j][i] ;
            for ( i=0; i<na; i++ )  Binv[j][i] -= top * Binv[k][i];
         }

      Bidx[k] = ell;
      // To save Bidx, x, Binv
      if ( ell >= TstPt && !RelTab[ell][FixPt] )
      {
         RelTab[ell][FixPt] = RelTab[FixPt][ell] = true;
         J[0] = FixPt;
         J[1] = ell; // J[0] < J[1] is assumed.
         L0.Add(J, na, Bidx, x, Binv);  // add 2 points
      }
   } // LP loop
} // EOF

////////////////////////////////////////

#include <string.h>
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "IT_LP.h"
#include "L0_IT_LP.h"

void dlp1_1pt_i
   (
   int ma, int na, double** a, int &bIdx, 
   double* c, int* LPidx, int &FixPt, int &TstPt, 
   int* Bidx, double* x, double** Binv, bool** RelTab 
   )
{
   bool   ibrk;
   int    i, i1, j, k, ell;
   double eps=1.0e-6, vmax, top, bot, sigj, dtmp;
 
   while(1) // LP loop
   {
      // STEP 1: FIND WHICH CONSTRAINT BEING KICKED OUT
      vmax = -DBL_MAX;
      for ( i=0; i<na; i++ ) 
      {
         dtmp = c[0] * Binv[i][0];
         for ( j=1; j<na; j++ )  dtmp += c[j] * Binv[i][j];
         if( dtmp > vmax )  { vmax = dtmp;  k = i; }
      }
 
      if (vmax < eps) return;      // out with optimal solution

      // STEP 2: FIND WHICH CONSTRAINT COMING IN
      sigj = DBL_MAX;
      ell = -1;
      for ( i1=0; i1<ma; i1++ )
      {
         i = LPidx[i1];
         ibrk = false;
         for ( j=0; j<na; j++ )  
            if ( Bidx[j] == i )  { ibrk = true;  break; }
         if ( ibrk )  continue;

         bot = a[i][0]*Binv[k][0];
         for ( j=1; j<na; j++ )  bot += a[i][j]*Binv[k][j];
         if ( bot >= -eps )  continue;

         top = -a[i][bIdx];
         for ( j=0; j<na; j++ )  top += a[i][j]*x[j];
         top /= bot;
         if ( top >= sigj )  continue;

         sigj = top;
         ell = i;
      }

      if ( ell < 0 )
      {
         cerr << "dlp1_1pt_i.cpp: LP unbounded\n\a";
         abort();       
      }

      // STEP 3: UPDATE x, Bidx, Binv
      for ( i=0; i<na; i++ )  x[i] -= sigj * Binv[k][i];

      top = a[ell][0] * Binv[k][0];
      for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[k][i];
      if( fabs(top) < eps )
      {
         cerr << "dlp1_1pt_i.cpp: Base matrix is singular\a\n";
         abort();
      }
      top = 1.0 / top;
      for ( i=0; i<na; i++ )  Binv[k][i] *= top;
      for ( j=0; j<na; j++ )
         if ( j != k )
         {
            top = a[ell][0] * Binv[j][0];
            for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[j][i];
            for ( i=0; i<na; i++ )  Binv[j][i] -= top * Binv[k][i];
         }

      Bidx[k] = ell;

      if( ell >= TstPt ) 
      {
         for ( i=0; i<na; i++ )
            if ( i != k && !RelTab[ell][Bidx[i]] )
               RelTab[ell][Bidx[i]] = RelTab[Bidx[i]][ell] = true;

         if ( !RelTab[ell][FixPt] ) 
            RelTab[ell][FixPt] = RelTab[FixPt][ell] = true;
      }
   } // LP loop
} // EOF

////////////////////////////////////////

#include <string.h>
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "IT_LP.h"
#include "L0_IT_LP.h"

void dlp1_1pt_s
   (
   int ma, int na, double** a, int &bIdx, 
   double* c, int* LPidx, int &FixPt, int &TstPt, 
   int* Bidx, double* x, double** Binv, bool** RelTab,
   L0_IML &L0
   )
{
   bool   ibrk;
   int    i, i1, j, k, ell;
   double eps=1.0e-6, vmax, top, bot, sigj, dtmp ;
 
   extern int* J;  //   int* J = new int [2];

   while(1) // LP loop
   {
      // STEP 1: FIND WHICH CONSTRAINT BEING KICKED OUT
      vmax = -DBL_MAX;
      for ( i=0; i<na; i++ ) 
      {
         dtmp = c[0]*Binv[i][0];
         for ( j=1; j<na; j++ )  dtmp += c[j]*Binv[i][j];
         if ( dtmp > vmax )  { vmax = dtmp;  k = i; }
      }
 
      if ( vmax < eps ) return;  // out with optimal solution

      // STEP 2: FIND WHICH CONSTRAINT COMING IN
      sigj = DBL_MAX;
      ell = -1;
      for ( i1=0; i1<ma; i1++ )
      {
         i = LPidx[i1];
         ibrk = false;
         for ( j=0; j<na; j++ )
            if ( Bidx[j] == i ) { ibrk = true;  break; }
         if ( ibrk )  continue;

         bot = a[i][0]*Binv[k][0];
         for ( j=1; j<na; j++ )  bot += a[i][j]*Binv[k][j];
         if ( bot >= -eps )  continue;

         top = -a[i][bIdx];
         for ( j=0; j<na; j++ )  top += a[i][j]*x[j];
         top /= bot;
         if ( top >= sigj )  continue;

         sigj = top;
         ell = i;
      }

      if ( ell < 0 )
      {
         cerr << "dlp1_1pt_s.cpp: LP unbounded\n\a";
         abort();       
      }

      // STEP 3: UPDATE x, Bidx, Binv
      for ( i=0; i<na; i++ )  x[i] -= sigj * Binv[k][i];

      top = a[ell][0] * Binv[k][0];
      for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[k][i];
      if ( fabs(top) < eps )
      {
         cerr << "dlp1_1pt_s.cpp: Base matrix is singular\a\n";
         abort();
      }
      top = 1.0 / top;
      for ( i=0; i<na; i++ )  Binv[k][i] *= top;
      for ( j=0; j<na; j++ )
         if ( j != k )
         {
            top = a[ell][0] * Binv[j][0];
            for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[j][i];
            for ( i=0; i<na; i++ )  Binv[j][i] -= top * Binv[k][i];
         }

      Bidx[k] = ell;

      // To save Bidx, x, Binv
      if ( ell >= TstPt && !RelTab[ell][FixPt] )
      {
         RelTab[ell][FixPt] = RelTab[FixPt][ell] = true;
         J[0] = FixPt;
         J[1] = ell;    // J[0] < J[1] is assumed.
         L0.Add(J, na, Bidx, x, Binv);  // add 2 points
      }
   } // LP loop
} // EOF

#endif
