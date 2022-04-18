#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* min/test.h
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

gsl_function create_function (MpIeee(*f)(MpIeee, void *));

void
test_f_e (const gsl_min_fminimizer_type * T, 
          const char * description, gsl_function *f,
          MpIeee lower_bound, MpIeee minimum, MpIeee upper_bound, 
          MpIeee correct_minimum);

void
test_f (const gsl_min_fminimizer_type * T, 
        const char * description, gsl_function *f,
        MpIeee lower_bound, MpIeee middle, MpIeee upper_bound, 
        MpIeee correct_minimum);

int
 test_bracket(const char * description,gsl_function *f,MpIeee lower_bound, 
              MpIeee upper_bound, unsigned int  max);

MpIeee f_cos(MpIeee x, void * p);
MpIeee func1(MpIeee x, void * p);
MpIeee func2(MpIeee x, void * p);
MpIeee func3(MpIeee x, void * p);
MpIeee func4(MpIeee x, void * p);
