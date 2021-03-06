#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_olver.h
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

#ifndef BESSEL_OLVER_H_
#define BESSEL_OLVER_H_

#include <gsl/gsl_sf_result.h>

int  gsl_sf_bessel_Jnu_asymp_Olver_e(MpIeee nu, MpIeee x, gsl_sf_result * result);
int  gsl_sf_bessel_Ynu_asymp_Olver_e(MpIeee nu, MpIeee x, gsl_sf_result * result);

MpIeee gsl_sf_bessel_Olver_zofmzeta(MpIeee minus_zeta);


#endif /* !BESSEL_OLVER_H_ */
