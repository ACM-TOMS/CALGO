#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_temme.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman */

#ifndef BESSEL_TEMME_H_
#define BESSEL_TEMME_H_

#include <gsl/gsl_sf_result.h>


int
 gsl_sf_bessel_Y_temme(const MpIeee nu, const MpIeee x,
                      gsl_sf_result * Y_nu,
                      gsl_sf_result * Y_nup1);

int
 gsl_sf_bessel_K_scaled_temme(const MpIeee nu, const MpIeee x,
                             MpIeee * K_nu, MpIeee * K_nup1, MpIeee * Kp_nu);


#endif /* !BESSEL_TEMME_H_ */
