#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/covar_source.c
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

static MpIeee FUNCTION(compute,covariance) (const BASE data1[], const size_t stride1,
                              const BASE data2[], const size_t stride2,
                              const size_t n, 
                              const MpIeee mean1, const MpIeee mean2);

static MpIeee FUNCTION(compute,covariance) (const BASE data1[], const size_t stride1,
                              const BASE data2[], const size_t stride2,
                              const size_t n, 
                              const MpIeee mean1, const MpIeee mean2)
{
  /* takes a dataset and finds the covariance */

  MpIeee covariance=  MpIeee( "0" ) ;

  size_t i;

  /* find the sum of the squares */
  for (i = 0; i < n; i++)
    {
      const MpIeee delta1=  (data1[i * stride1] - mean1);
      const MpIeee delta2=  (data2[i * stride2] - mean2);
      covariance += (delta1 * delta2 - covariance) / (i + MpIeee( "1" ));
    }

  return covariance ;
}

MpIeee FUNCTION(gsl_stats,covariance_m) (const BASE data1[], const size_t stride1, 
                                  const BASE data2[], const size_t stride2, 
                                  const size_t n, 
                                  const MpIeee mean1, const MpIeee mean2)
{
  const MpIeee covariance=  FUNCTION(compute,covariance) (data1, stride1,
                                                          data2, stride2,
                                                          n, 
                                                          mean1, mean2);
  
  return covariance * ((MpIeee)n / (MpIeee)(n - MpIeee( "1" )));
}

MpIeee FUNCTION(gsl_stats,covariance) (const BASE data1[], const size_t stride1,
                                const BASE data2[], const size_t stride2,
                                const size_t n)
{
  const MpIeee mean1=  FUNCTION(gsl_stats,mean) (data1, stride1, n);
  const MpIeee mean2=  FUNCTION(gsl_stats,mean) (data2, stride2, n);

  return FUNCTION(gsl_stats,covariance_m)(data1, stride1, 
                                          data2, stride2, 
                                          n, 
                                          mean1, mean2);
}


