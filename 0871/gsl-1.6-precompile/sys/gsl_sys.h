#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* sys/gsl_sys.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
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

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const MpIeee x);
MpIeee gsl_expm1(const MpIeee x);
MpIeee gsl_hypot(const MpIeee x, const MpIeee y);
MpIeee gsl_acosh(const MpIeee x);
MpIeee gsl_asinh(const MpIeee x);
MpIeee gsl_atanh(const MpIeee x);

int  gsl_isnan(const MpIeee x);
int  gsl_isinf(const MpIeee x);
int  gsl_finite(const MpIeee x);

MpIeee gsl_nan(void);
MpIeee gsl_posinf(void);
MpIeee gsl_neginf(void);
MpIeee gsl_fdiv(const MpIeee x, const MpIeee y);

MpIeee gsl_coerce_double(const MpIeee x);
MpIeee gsl_coerce_float(const MpIeee x);
MpIeee gsl_coerce_long_double(const MpIeee x);

MpIeee gsl_ldexp(const MpIeee x, const int e);
MpIeee gsl_frexp(const MpIeee x, int  * e);

int  gsl_fcmp(const MpIeee x1, const MpIeee x2, const MpIeee epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */
