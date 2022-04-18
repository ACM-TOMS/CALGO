#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* fit/gsl_fit.h
 * 
 * Copyright (C) 2000 Brian Gough
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

#ifndef __GSL_FIT_H__
#define __GSL_FIT_H__

#include <stdlib.h>
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

int gsl_fit_linear (const MpIeee * x, const size_t xstride,
                    const MpIeee * y, const size_t ystride,
                    const size_t n,
                    MpIeee * c0, MpIeee * c1, 
                    MpIeee * cov00, MpIeee * cov01, MpIeee * cov11, 
                    MpIeee * sumsq);


int  gsl_fit_wlinear(const MpIeee * x, const size_t xstride,
                     const MpIeee * w, const size_t wstride,
                     const MpIeee * y, const size_t ystride,
                     const size_t n,
                     MpIeee * c0, MpIeee * c1, 
                     MpIeee * cov00, MpIeee * cov01, MpIeee * cov11, 
                     MpIeee * chisq);

int
 gsl_fit_linear_est(const MpIeee x, 
                    const MpIeee c0, const MpIeee c1, 
                    const MpIeee c00, const MpIeee c01, const MpIeee c11,
                    MpIeee *y, MpIeee *y_err);


int  gsl_fit_mul(const MpIeee * x, const size_t xstride,
                 const MpIeee * y, const size_t ystride,
                 const size_t n,
                 MpIeee * c1, 
                 MpIeee * cov11, 
                 MpIeee * sumsq);

int  gsl_fit_wmul(const MpIeee * x, const size_t xstride,
                  const MpIeee * w, const size_t wstride,
                  const MpIeee * y, const size_t ystride,
                  const size_t n,
                  MpIeee * c1, 
                  MpIeee * cov11, 
                  MpIeee * sumsq);


int
 gsl_fit_mul_est(const MpIeee x, 
                 const MpIeee c1, 
                 const MpIeee c11,
                 MpIeee *y, MpIeee *y_err);

__END_DECLS

#endif /* __GSL_FIT_H__ */
