#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/wvariance_source.c
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

static MpIeee FUNCTION(compute,wvariance) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const MpIeee wmean);

static MpIeee FUNCTION(compute,factor) (const BASE w[], const size_t wstride, const size_t n);

static MpIeee FUNCTION(compute,wvariance) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const MpIeee wmean)
{
  /* takes a dataset and finds the weighted variance */

  MpIeee wvariance=  MpIeee( "0" ) ;
  MpIeee W=  MpIeee( "0" );

  size_t i;

  /* find the sum of the squares */
  for (i = 0; i < n; i++)
    {
      MpIeee wi=  w[i * wstride];

      if (wi > MpIeee( "0" )) {
        const MpIeee delta=  (data[i * stride] - wmean);
        W += wi ;
        wvariance += (delta * delta - wvariance) * (wi / W);
      }
    }

  return wvariance ;
}

static MpIeee FUNCTION(compute,factor) (const BASE w[], const size_t wstride, const size_t n)
{
  /* Find the factor ``N/(N-1)'' which multiplies the raw std dev */

  MpIeee a=  MpIeee( "0" ) ;
  MpIeee b=  MpIeee( "0" );
  MpIeee factor;

  size_t i;

  /* find the sum of the squares */
  for (i = 0; i < n; i++)
    {
      MpIeee wi=  w[i * wstride];

      if (wi > MpIeee( "0" ))
        {
          a += wi ;
          b += wi * wi ;
        }
    }

  factor = (a*a) / ((a*a) - b);

  return factor ;
}

MpIeee FUNCTION(gsl_stats,wvariance_with_fixed_mean) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const MpIeee wmean)
{
  const MpIeee wvariance=  FUNCTION(compute,wvariance) (w, wstride, data, stride, n, wmean);
  return wvariance;
}

MpIeee FUNCTION(gsl_stats,wsd_with_fixed_mean) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const MpIeee wmean)
{
  const MpIeee wvariance=  FUNCTION(compute,wvariance) (w, wstride, data, stride, n, wmean);
  const MpIeee wsd=  sqrt (wvariance);

  return wsd;
}


MpIeee FUNCTION(gsl_stats,wvariance_m) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const MpIeee wmean)
{
  const MpIeee variance=  FUNCTION(compute,wvariance) (w, wstride, data, stride, n, wmean);
  const MpIeee scale=  FUNCTION(compute,factor)(w, wstride, n);
  
  return scale * variance;
}

MpIeee FUNCTION(gsl_stats,wsd_m) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n, const MpIeee wmean)
{
  const MpIeee variance=  FUNCTION(compute,wvariance) (w, wstride, data, stride, n, wmean);
  const MpIeee scale=  FUNCTION(compute,factor)(w, wstride, n);
  const MpIeee wsd=  sqrt(scale * variance) ;
  
  return wsd;
}

MpIeee FUNCTION(gsl_stats,wsd) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n)
{
  const MpIeee wmean=  FUNCTION(gsl_stats,wmean) (w, wstride, data, stride, n);
  return FUNCTION(gsl_stats,wsd_m) (w, wstride, data, stride, n, wmean) ;
}

MpIeee FUNCTION(gsl_stats,wvariance) (const BASE w[], const size_t wstride, const BASE data[], const size_t stride, const size_t n)
{
  const MpIeee wmean=  FUNCTION(gsl_stats,wmean) (w, wstride, data, stride, n);
  return FUNCTION(gsl_stats,wvariance_m)(w, wstride, data, stride, n, wmean);
}
