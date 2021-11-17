#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/variance_source.c
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

static MpIeee FUNCTION(compute,variance) (const BASE data[], const size_t stride, const size_t n, const MpIeee mean);

static MpIeee FUNCTION(compute,variance) (const BASE data[], const size_t stride, const size_t n, const MpIeee mean)
{
  /* takes a dataset and finds the variance */

  MpIeee variance=  MpIeee( "0" ) ;

  size_t i;

  /* find the sum of the squares */
  for (i = 0; i < n; i++)
    {
      const MpIeee delta=  (data[i * stride] - mean);
      variance += (delta * delta - variance) / (i + MpIeee( "1" ));
    }

  return variance ;
}


MpIeee FUNCTION(gsl_stats,variance_with_fixed_mean) (const BASE data[], const size_t stride, const size_t n, const MpIeee mean)
{
  const MpIeee variance=  FUNCTION(compute,variance) (data, stride, n, mean);
  return variance;
}

MpIeee FUNCTION(gsl_stats,sd_with_fixed_mean) (const BASE data[], const size_t stride, const size_t n, const MpIeee mean)
{
  const MpIeee variance=  FUNCTION(compute,variance) (data, stride, n, mean);
  const MpIeee sd=  sqrt (variance);

  return sd;
}



MpIeee FUNCTION(gsl_stats,variance_m) (const BASE data[], const size_t stride, const size_t n, const MpIeee mean)
{
  const MpIeee variance=  FUNCTION(compute,variance) (data, stride, n, mean);
  
  return variance * ((MpIeee)n / (MpIeee)(n - MpIeee( "1" )));
}

MpIeee FUNCTION(gsl_stats,sd_m) (const BASE data[], const size_t stride, const size_t n, const MpIeee mean)
{
  const MpIeee variance=  FUNCTION(compute,variance) (data, stride, n, mean);
  const MpIeee sd=  sqrt (variance * ((const MpIeee)n / (const MpIeee)(n - 1)));

  return sd;
}

MpIeee FUNCTION(gsl_stats,variance) (const BASE data[], const size_t stride, const size_t n)
{
  const MpIeee mean=  FUNCTION(gsl_stats,mean) (data, stride, n);
  return FUNCTION(gsl_stats,variance_m)(data, stride, n, mean);
}

MpIeee FUNCTION(gsl_stats,sd) (const BASE data[], const size_t stride, const size_t n)
{
  const MpIeee mean=  FUNCTION(gsl_stats,mean) (data, stride, n);
  return FUNCTION(gsl_stats,sd_m) (data, stride, n, mean);
}
