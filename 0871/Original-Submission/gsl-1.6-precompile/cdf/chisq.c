#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* cdf/cdf_chisq.c
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA.
 */

#include <config.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>

MpIeee gsl_cdf_chisq_P(const MpIeee x, const MpIeee nu)
{
  return gsl_cdf_gamma_P (x, nu / MpIeee( "2" ), MpIeee( "2.0" ));
}

MpIeee gsl_cdf_chisq_Q(const MpIeee x, const MpIeee nu)
{
  return gsl_cdf_gamma_Q (x, nu / MpIeee( "2" ), MpIeee( "2.0" ));
}
