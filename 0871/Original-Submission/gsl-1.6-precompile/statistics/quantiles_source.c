#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/quantiles_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Jim Davies, Brian Gough
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


MpIeee FUNCTION(gsl_stats,quantile_from_sorted_data) (const BASE sorted_data[], 
                                               const size_t stride,
                                               const size_t n,
                                               const MpIeee f)
{
  const MpIeee index=  f * (n - 1) ;
  const size_t lhs = (int)index.toInt() ;
  const MpIeee delta=  index - lhs ;
  MpIeee result;

  if (n == 0)
    return MpIeee( "0.0" ) ;

  if (lhs == n - 1)
    {
      result = sorted_data[lhs * stride] ;
    }
  else 
    {
      result = (MpIeee( "1" ) - delta) * sorted_data[lhs * stride] + delta * sorted_data[(lhs + 1) * stride] ;
    }

  return result ;
}
