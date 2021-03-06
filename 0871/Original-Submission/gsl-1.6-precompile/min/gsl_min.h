#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* min/gsl_min.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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

#ifndef __GSL_MIN_H__
#define __GSL_MIN_H__

#include <stdlib.h>
#include <gsl/gsl_types.h>
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

typedef struct
  {
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_function * f, MpIeee x_minimum, MpIeee f_minimum, MpIeee x_lower, MpIeee f_lower, MpIeee x_upper, MpIeee f_upper);
    int (*iterate) (void *state, gsl_function * f, MpIeee * x_minimum, MpIeee * f_minimum, MpIeee * x_lower, MpIeee * f_lower, MpIeee * x_upper, MpIeee * f_upper);
  }
gsl_min_fminimizer_type;

typedef struct
  {
    const gsl_min_fminimizer_type * type;
    gsl_function * function ;
    MpIeee x_minimum;
    MpIeee x_lower;
    MpIeee x_upper;
    MpIeee f_minimum;MpIeee  f_lower;MpIeee  f_upper;
    void *state;
  }
gsl_min_fminimizer;

gsl_min_fminimizer *
gsl_min_fminimizer_alloc (const gsl_min_fminimizer_type * T) ;
                                      
void gsl_min_fminimizer_free (gsl_min_fminimizer * s);

int  gsl_min_fminimizer_set(gsl_min_fminimizer * s, 
                            gsl_function * f, MpIeee x_minimum, 
                            MpIeee x_lower, MpIeee x_upper);

int  gsl_min_fminimizer_set_with_values(gsl_min_fminimizer * s, 
                                        gsl_function * f, 
                                        MpIeee x_minimum, MpIeee f_minimum,
                                        MpIeee x_lower, MpIeee f_lower,
                                        MpIeee x_upper, MpIeee f_upper);

int  gsl_min_fminimizer_iterate(gsl_min_fminimizer * s);

const char * gsl_min_fminimizer_name (const gsl_min_fminimizer * s);

MpIeee gsl_min_fminimizer_x_minimum(const gsl_min_fminimizer * s);
MpIeee gsl_min_fminimizer_x_lower(const gsl_min_fminimizer * s);
MpIeee gsl_min_fminimizer_x_upper(const gsl_min_fminimizer * s);
MpIeee gsl_min_fminimizer_f_minimum(const gsl_min_fminimizer * s);
MpIeee gsl_min_fminimizer_f_lower(const gsl_min_fminimizer * s);
MpIeee gsl_min_fminimizer_f_upper(const gsl_min_fminimizer * s);

/* Deprecated, use x_minimum instead */
MpIeee gsl_min_fminimizer_minimum(const gsl_min_fminimizer * s);

int
 gsl_min_test_interval(MpIeee x_lower, MpIeee x_upper, MpIeee epsabs, MpIeee epsrel);

GSL_VAR const gsl_min_fminimizer_type  * gsl_min_fminimizer_goldensection;
GSL_VAR const gsl_min_fminimizer_type  * gsl_min_fminimizer_brent;

typedef
int (*gsl_min_bracketing_function)(gsl_function *f,
                                   MpIeee *x_minimum,MpIeee * f_minimum,
                                   MpIeee *x_lower, MpIeee * f_lower,
                                   MpIeee *x_upper, MpIeee * f_upper,
                                   size_t eval_max);

int 
 gsl_min_find_bracket(gsl_function *f,MpIeee *x_minimum,MpIeee * f_minimum,
                     MpIeee *x_lower, MpIeee * f_lower,
                     MpIeee *x_upper, MpIeee * f_upper,
                     size_t eval_max);

__END_DECLS

#endif /* __GSL_MIN_H__ */
