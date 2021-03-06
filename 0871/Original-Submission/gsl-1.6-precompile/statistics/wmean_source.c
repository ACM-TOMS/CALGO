#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/wmean_source.c
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

MpIeee FUNCTION(gsl_stats, wmean) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t size)
{
  /* Compute the weighted arithmetic mean M of a dataset using the
     recurrence relation

     M(n) = M(n-1) + (data[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
     W(n) = W(n-1) + w(n)

   */

  MpIeee wmean=  MpIeee( "0" );
  MpIeee W=  MpIeee( "0" );

  size_t i;

  for (i = 0; i < size; i++)
    {
      MpIeee wi=  w[i * wstride];

      if (wi > MpIeee( "0" ))
        {
          W += wi;
          wmean += (data[i * stride] - wmean) * (wi / W);
        }
    }

  return wmean;
}
