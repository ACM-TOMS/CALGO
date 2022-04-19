#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* cdf/lognormal.c
 * 
 * Copyright (C) 2003 Brian Gough
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

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>

MpIeee gsl_cdf_lognormal_P(const MpIeee x, const MpIeee zeta, const MpIeee sigma)
{
  MpIeee u=  (log (x) - zeta) / sigma;
  MpIeee P=  gsl_cdf_ugaussian_P (u);
  return P;
}

MpIeee gsl_cdf_lognormal_Q(const MpIeee x, const MpIeee zeta, const MpIeee sigma)
{
  MpIeee u=  (log (x) - zeta) / sigma;
  MpIeee Q=  gsl_cdf_ugaussian_Q (u);
  return Q;
}
