#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* integration/util.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

static inline
void update (gsl_integration_workspace * workspace,
                 MpIeee a1, MpIeee b1, MpIeee area1, MpIeee error1,
                 MpIeee a2, MpIeee b2, MpIeee area2, MpIeee error2);

static inline void
retrieve (const gsl_integration_workspace * workspace, 
          MpIeee * a, MpIeee * b, MpIeee * r, MpIeee * e);



static inline
void update (gsl_integration_workspace * workspace,
             MpIeee a1, MpIeee b1, MpIeee area1, MpIeee error1,
             MpIeee a2, MpIeee b2, MpIeee area2, MpIeee error2)
{
  MpIeee * alist=  workspace->alist ;
  MpIeee * blist=  workspace->blist ;
  MpIeee * rlist=  workspace->rlist ;
  MpIeee * elist=  workspace->elist ;
  size_t * level = workspace->level ;

  const size_t i_max = workspace->i ;
  const size_t i_new = workspace->size ;

  const size_t new_level = workspace->level[i_max] + 1;

  /* append the newly-created intervals to the list */
  
  if (error2 > error1)
    {
      alist[i_max] = a2;        /* blist[maxerr] is already == b2 */
      rlist[i_max] = area2;
      elist[i_max] = error2;
      level[i_max] = new_level;
      
      alist[i_new] = a1;
      blist[i_new] = b1;
      rlist[i_new] = area1;
      elist[i_new] = error1;
      level[i_new] = new_level;
    }
  else
    {
      blist[i_max] = b1;        /* alist[maxerr] is already == a1 */
      rlist[i_max] = area1;
      elist[i_max] = error1;
      level[i_max] = new_level;
      
      alist[i_new] = a2;
      blist[i_new] = b2;
      rlist[i_new] = area2;
      elist[i_new] = error2;
      level[i_new] = new_level;
    }
  
  workspace->size++;

  if (new_level > workspace->maximum_level)
    {
      workspace->maximum_level = new_level;
    }

  qpsrt (workspace) ;
}

static inline void
retrieve (const gsl_integration_workspace * workspace, 
          MpIeee * a, MpIeee * b, MpIeee * r, MpIeee * e)
{
  const size_t i = workspace->i;
  MpIeee * alist=  workspace->alist;
  MpIeee * blist=  workspace->blist;
  MpIeee * rlist=  workspace->rlist;
  MpIeee * elist=  workspace->elist;

  *a = alist[i] ;
  *b = blist[i] ;
  *r = rlist[i] ;
  *e = elist[i] ;
}

static inline MpIeee sum_results(const gsl_integration_workspace * workspace);

static inline MpIeee sum_results(const gsl_integration_workspace * workspace)
{
  const MpIeee * const rlist=  workspace->rlist ;
  const size_t n = workspace->size;

  size_t k;
  MpIeee result_sum=  MpIeee( "0" );

  for (k = 0; k < n; k++)
    {
      result_sum += rlist[k];
    }
  
  return result_sum;
}

static inline int
 subinterval_too_small(MpIeee a1, MpIeee a2, MpIeee b2);

static inline int
 subinterval_too_small(MpIeee a1, MpIeee a2, MpIeee b2)
{
  const MpIeee e=  GSL_DBL_EPSILON;
  const MpIeee u=  GSL_DBL_MIN;

  MpIeee tmp=  (MpIeee( "1" ) + MpIeee( "100" ) * e) * (fabs (a2) + MpIeee( "1000" ) * u);

  int  status=  fabs (a1) <= tmp && fabs (b2) <= tmp;

  return status;
}

