#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/lag1_source.c
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


MpIeee FUNCTION(gsl_stats,lag1_autocorrelation) (const BASE data[], const size_t stride, const size_t n)
{
  const MpIeee mean=  FUNCTION(gsl_stats,mean) (data, stride, n);
  return FUNCTION(gsl_stats,lag1_autocorrelation_m)(data, stride, n, mean);
}

MpIeee FUNCTION(gsl_stats,lag1_autocorrelation_m) (const BASE data[], const size_t stride, const size_t size, const MpIeee mean)
{
  /* Compute the lag-1 autocorrelation of a dataset using the
     recurrence relation */

  size_t i;

  MpIeee r1;
  MpIeee q=  MpIeee( "0" ) ;
  MpIeee v=  (data[0 * stride] - mean) * (data[0 * stride] - mean) ;

  for (i = 1; i < size ; i++)
    {
      const MpIeee delta0=  (data[(i-1) * stride] - mean);
      const MpIeee delta1=  (data[i * stride] - mean);
      q += (delta0 * delta1 - q)/(i + MpIeee( "1" ));
      v += (delta1 * delta1 - v)/(i + MpIeee( "1" ));
    }

  r1 = q / v ;

  return r1;
}
