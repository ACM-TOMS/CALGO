#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* deriv/gsl_deriv.h
 * 
 * Copyright (C) 2000 David Morrison
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

#ifndef __GSL_DERIV_H__
#define __GSL_DERIV_H__
#include <gsl/gsl_math.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

int gsl_deriv_central (const gsl_function *f,
                       MpIeee x, MpIeee h,
                       MpIeee *result, MpIeee *abserr);

int  gsl_deriv_backward(const gsl_function *f,
                        MpIeee x, MpIeee h,
                        MpIeee *result, MpIeee *abserr);

int  gsl_deriv_forward(const gsl_function *f,
                       MpIeee x, MpIeee h,
                       MpIeee *result, MpIeee *abserr);

__END_DECLS

#endif /* __GSL_DERIV_H__ */
