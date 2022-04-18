#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/wskew_source.c
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

MpIeee FUNCTION(gsl_stats,wskew) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n)
{
  const MpIeee wmean=  FUNCTION(gsl_stats,wmean)(w, wstride, data, stride, n);
  const MpIeee wsd=  FUNCTION(gsl_stats,wsd_m)(w, wstride, data, stride, n, wmean);
  return FUNCTION(gsl_stats,wskew_m_sd)(w, wstride, data, stride, n, wmean, wsd);
}
    
MpIeee FUNCTION(gsl_stats,wskew_m_sd) (const BASE w[], const size_t wstride, 
                                const BASE data[], 
                                const size_t stride, const size_t n,
                                const MpIeee wmean, const MpIeee wsd)
{
  /* Compute the weighted skewness of a dataset */

  MpIeee wskew=  MpIeee( "0" );
  MpIeee W=  MpIeee( "0" );

  size_t i;

  /* find the sum of the cubed deviations, normalized by the sd. */

  /* we use a recurrence relation to stably update a running value so
     there aren't any large sums that can overflow */

  for (i = 0; i < n; i++)
    {
      MpIeee wi=  w[i * wstride];
      
      if (wi > MpIeee( "0" )) {
        const MpIeee x=  (data[i * stride] - wmean) / wsd;
        W += wi ;
        wskew += (x * x * x - wskew) * (wi / W);
      }
    }

  return wskew;
}

