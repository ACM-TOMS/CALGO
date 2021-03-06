#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* fft/compare.h
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

int
 FUNCTION(compare_complex,results) (const char *name_a, const BASE a[],
                                   const char *name_b, const BASE b[],
                                   size_t stride, size_t n, 
                                   const MpIeee allowed_ticks);

int
 FUNCTION(compare_real,results) (const char *name_a, const BASE a[],
                                   const char *name_b, const BASE b[],
                                   size_t stride, size_t n, 
                                   const MpIeee allowed_ticks);
