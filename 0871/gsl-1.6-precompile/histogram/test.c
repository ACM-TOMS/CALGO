#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* histogram/test.c
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

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

void test1d (void);
void test2d (void);
void test1d_resample (void);
void test2d_resample (void);
void test1d_trap (void);
void test2d_trap (void);

int
main (void)
{
  test1d();
  test2d();
  test1d_resample();
  test2d_resample();
  test1d_trap();
  test2d_trap();
  
  exit (gsl_test_summary ());
}
