#ifndef _Extention_LPCODE_
#define _Extention_LPCODE_  

#include <string.h>
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "IT_LP.h"

void dnulp2_a 
   (
   int ma, int na, double** a, int &bIdx, double* c,
   int* Bidx, double* x, double** Binv, 
   int &info
   )
{

   bool   ibrk;
   int    ell, k, k_, i, j, flag;
   double eps=1.0e-6, sigj, vmax, bot, top, dtmp;

   extern double* v;   // double* v = new double [na];   

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
      for ( i=0; i<ma; i++ )
      {
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
         for ( i=0; i<ma; i++ )
         {
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

void dlp2_1pts
   (
   int ma, int &na, double** a, int &bIdx,
   double* c, int &TstPt,
   int* Bidx, double* x, double** Binv,
   bool* PtIn, IT_LP &ItLp
   )

{
   bool   ibrk;
   int    ell, k, i, j;
   double eps=1.0e-6, sigj, vmax, bot, top, dtmp;

   while(1) // LP loop
   {
      // STEP 1: FIND WHICH CONSTRAINT BEING KICKED OUT
      vmax = -DBL_MAX;
      for ( i=0; i<na; i++ )
         if ( Bidx[i] != -1 )
         {
            dtmp = c[0]*Binv[i][0];
            for ( j=1; j<na; j++ )  dtmp += c[j]*Binv[i][j];
            if( dtmp > vmax )  { vmax = dtmp;  k = i; }
         }

      if (vmax < eps) return;      // out with optimal solution

      // STEP 2: FIND WHICH CONSTRAINT COMING IN
      sigj = DBL_MAX;
      ell = -1;
      for ( i=0; i<ma; i++ )
      {
         ibrk = false;
         for ( j=0; j<na; j++ )
            if ( Bidx[j] == i ) { ibrk = true;  break; }
         if (ibrk)  continue;

         bot = a[i][0] * Binv[k][0];
         for ( j=1; j<na; j++ ) bot += a[i][j] * Binv[k][j];
         if ( bot >= -eps )  continue;

         top = -a[i][bIdx];
         for ( j=0; j<na; j++ )  top += a[i][j] * x[j];
         top /= bot;
         if ( top >= sigj )  continue;

         sigj = top;
         ell = i;
      }

      if ( ell < 0 )
      {
         cerr << "dlp2_1pts.cpp: LP unbounded\n\a";
         abort();
      }

      // STEP 3: UPDATE x, Bidx, Binv
      for ( i=0; i<na; i++ )  x[i] -= sigj * Binv[k][i];

      top = a[ell][0] * Binv[k][0];
      for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[k][i];
      if( fabs(top) < eps )
      {
         cerr << "dlp2_1pts.cpp: Base matrix is singular\a\n";
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
      if ( ell >= TstPt && !PtIn[ell] )
      {
         PtIn[ell] = true;
         ItLp.Add(ell, na, Bidx, x, Binv);
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

void dlp1_1pts
   (
   int ma, int &na, double** a, int &bIdx,
   double* c, int &TstPt,
   int* Bidx, double* x, double** Binv,
   bool* PtIn, IT_LP &ItLp
   )
{
   int    i, j, k, ell;
   double eps=1.0e-6, vmax, top, bot, sigj, dtmp;
   bool   ibrk;

   while(1) // LP loop
   {
      // STEP 1: FIND WHICH CONSTRAINT BEING KICKED OUT
      vmax = -DBL_MAX;
      for(i=0; i<na; i++)
      {
         dtmp = c[0]*Binv[i][0];
         for ( j=1; j<na; j++ )  dtmp += c[j]*Binv[i][j];
         if(dtmp > vmax)  { vmax = dtmp;  k = i; }
      }

      if (vmax < eps) return;      // out with optimal solution

      // STEP 2: FIND WHICH CONSTRAINT COMING IN
      sigj = DBL_MAX;
      ell = -1;
      for ( i=0; i<ma; i++ )
      {
         ibrk = false;
         for( j=0; j<na; j++ )
            if (Bidx[j] == i) {ibrk = true;  break;}
         if (ibrk)  continue;

         bot = a[i][0] * Binv[k][0];
         for ( j=1; j<na; j++ ) bot += a[i][j] * Binv[k][j];
         if ( bot > -eps )  continue;

         top = -a[i][bIdx];
         for ( j=0; j<na; j++ )  top += a[i][j] * x[j];
         top /= bot;
         if (top > sigj)  continue;

         sigj = top;
         ell = i;
      }
      if ( ell < 0 )
      {
         cerr << "dlp1_1pts.cpp: LP unbounded\n\a";
         abort();
      }

      // STEP 3: UPDATE x, Bidx, Binv
      for ( i=0; i<na; i++ )  x[i] -= sigj * Binv[k][i];

      top = a[ell][0] * Binv[k][0];
      for ( i=1; i<na; i++ )  top += a[ell][i] * Binv[k][i];

      if( fabs(top) < eps )
      {
         cerr << "dlp1_1pts.cpp: Base matrix is singular\a\n";
         abort();
      }
      top = 1.0 / top;
      for ( i=0; i<na; i++ )  Binv[k][i] *= top;
      for ( j=0; j<na; j++ )
      if ( j != k )
      {
         top = a[ell][0] * Binv[j][0];
         for ( i=1; i <na; i++ )  top += a[ell][i] * Binv[j][i] ;
         for ( i=0; i <na; i++ )  Binv[j][i] -= top * Binv[k][i];
      }

      Bidx[k] = ell;

      // To save Bidx, x, Binv

      if ( ell >= TstPt && !PtIn[ell] )
      {
         PtIn[ell] = true;
         ItLp.Add(ell, na, Bidx, x, Binv);
      }
   } // LP loop
} // EOF

#endif
