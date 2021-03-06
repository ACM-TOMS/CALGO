#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/ttest_source.c
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

MpIeee FUNCTION(gsl_stats,ttest) (const BASE data1[], 
                           const size_t stride1, const size_t n1, 
                           const BASE data2[],
                           const size_t stride2, const size_t n2)
{
  /* runs a t-test between two datasets representing independent
     samples. Tests to see if the difference between means of the
     samples is different from zero */

  /* find means for the two samples */
  const MpIeee mean1=  FUNCTION(gsl_stats,mean) (data1, stride1, n1);
  const MpIeee mean2=  FUNCTION(gsl_stats,mean) (data2, stride2, n2);

  /* find pooled variance for the two samples */
  const MpIeee pv=  FUNCTION(gsl_stats,pvariance) (data1, stride1, n1, data2, stride2, n2);

  /* calculate the t statistic */
  const MpIeee t=  (mean1 - mean2) / (sqrt (pv * ((1.0 / n1) + (1.0 / n2))));

  return t;
}

