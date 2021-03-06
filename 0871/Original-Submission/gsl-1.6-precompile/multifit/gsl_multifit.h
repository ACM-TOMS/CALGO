#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* multifit/gsl_multifit.h
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

#ifndef __GSL_MULTIFIT_H__
#define __GSL_MULTIFIT_H__

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

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

typedef struct 
{
  size_t n; /* number of observations */
  size_t p; /* number of parameters */
  gsl_matrix * A;
  gsl_matrix * Q;
  gsl_matrix * QSI;
  gsl_vector * S;
  gsl_vector * t;
  gsl_vector * xt;
  gsl_vector * D;
} 
gsl_multifit_linear_workspace;

gsl_multifit_linear_workspace *
gsl_multifit_linear_alloc (size_t n, size_t p);

void
gsl_multifit_linear_free (gsl_multifit_linear_workspace * work);

int
 gsl_multifit_linear(const gsl_matrix * X,
                     const gsl_vector * y,
                     gsl_vector * c,
                     gsl_matrix * cov,
                     MpIeee * chisq,
                     gsl_multifit_linear_workspace * work);

int
 gsl_multifit_wlinear(const gsl_matrix * X,
                      const gsl_vector * w,
                      const gsl_vector * y,
                      gsl_vector * c,
                      gsl_matrix * cov,
                      MpIeee * chisq,
                      gsl_multifit_linear_workspace * work);

__END_DECLS

#endif /* __GSL_MULTIFIT_H__ */
