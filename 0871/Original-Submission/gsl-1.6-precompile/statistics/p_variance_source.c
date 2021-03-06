#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/p_variance_source.c
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


MpIeee FUNCTION(gsl_stats,pvariance) (const BASE data1[], 
                               const size_t stride1, const size_t n1, 
                               const BASE data2[], 
                               const size_t stride2, const size_t n2)
{
  /* Find the pooled variance of two datasets */

  const MpIeee var1=  FUNCTION(gsl_stats,variance) (data1, stride1, n1);
  const MpIeee var2=  FUNCTION(gsl_stats,variance) (data2, stride2, n2);

  /* calculate the pooled variance */

  const MpIeee pooled_variance=  
    (((n1 - 1) * var1) + ((n2 - 1) * var2)) / (n1 + n2 - 2);

  return pooled_variance;
}

