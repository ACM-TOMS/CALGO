#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* roots/test.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Reid Priedhorsky, Brian Gough
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

gsl_function create_function (MpIeee(*f)(MpIeee, void *)) ;
gsl_function_fdf create_fdf (MpIeee(*f)(MpIeee, void *),
                             MpIeee(*df)(MpIeee, void *),
                             void (*fdf)(MpIeee, void *, MpIeee*, MpIeee*));

void
  test_macros (void);

void
  test_roots (void);

void
  test_poly (void);

void
test_f (const gsl_root_fsolver_type * T, 
        const char * description, gsl_function *f,
        MpIeee lower_bound, MpIeee upper_bound, MpIeee correct_root);

void
test_f_e (const gsl_root_fsolver_type * T, const char * description, 
          gsl_function *f,
          MpIeee lower_bound, MpIeee upper_bound, MpIeee correct_root);

void
test_fdf (const gsl_root_fdfsolver_type * T, const char * description, 
          gsl_function_fdf *fdf, MpIeee root, MpIeee correct_root);

void
test_fdf_e (const gsl_root_fdfsolver_type * T, const char * description, 
            gsl_function_fdf *fdf, MpIeee root, MpIeee correct_root);


void
  usage (void);

void
  error_handler (const char *reason, const char *file, int  line);

MpIeee func1(MpIeee x, void * p);

MpIeee func1_df(MpIeee x, void * p);

void
  func1_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime);

MpIeee func2(MpIeee x, void * p);

MpIeee func2_df(MpIeee x, void * p);

void
  func2_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime);

MpIeee func3(MpIeee x, void * p);

MpIeee func3_df(MpIeee x, void * p);

void
  func3_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime);

MpIeee func4(MpIeee x, void * p);

MpIeee func4_df(MpIeee x, void * p);

void
  func4_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime);

MpIeee func5(MpIeee x, void * p);

MpIeee func5_df(MpIeee x, void * p);

void
  func5_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime);

MpIeee func6(MpIeee x, void * p);

MpIeee func6_df(MpIeee x, void * p);

void
  func6_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime);

MpIeee sin_f(MpIeee x, void * p);

MpIeee sin_df(MpIeee x, void * p);

void
  sin_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime);

MpIeee cos_f(MpIeee x, void * p);

MpIeee cos_df(MpIeee x, void * p);

void
  cos_fdf (MpIeee x, void * p, MpIeee *y, MpIeee *yprime);
