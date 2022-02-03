#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* monte/gsl_monte_miser.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Michael Booth
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

/* Author: MJB */

#ifndef __GSL_MONTE_MISER_H__
#define __GSL_MONTE_MISER_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>

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

typedef struct {
  size_t min_calls;
  size_t min_calls_per_bisection;
  MpIeee dither;
  MpIeee estimate_frac;
  MpIeee alpha;
  size_t dim;
  int  estimate_style;
  int  depth;
  int  verbose;
  MpIeee * x;
  MpIeee * xmid;
  MpIeee * sigma_l;
  MpIeee * sigma_r;
  MpIeee * fmax_l;
  MpIeee * fmax_r;
  MpIeee * fmin_l;
  MpIeee * fmin_r;
  MpIeee * fsum_l;
  MpIeee * fsum_r;
  MpIeee * fsum2_l;
  MpIeee * fsum2_r;
  size_t * hits_l;
  size_t * hits_r;
} gsl_monte_miser_state; 

int  gsl_monte_miser_integrate(gsl_monte_function * f, 
                              const MpIeee xl[], const MpIeee xh[], 
                              size_t dim, size_t calls, 
                              gsl_rng *r, 
                              gsl_monte_miser_state* state,
                              MpIeee *result, MpIeee *abserr);

gsl_monte_miser_state* gsl_monte_miser_alloc(size_t dim);

int  gsl_monte_miser_init(gsl_monte_miser_state* state);

void gsl_monte_miser_free(gsl_monte_miser_state* state);


__END_DECLS

#endif /* __GSL_MONTE_MISER_H__ */
