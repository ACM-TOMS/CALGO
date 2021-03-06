#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* fft/hc_pass.h
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

#include "complex_internal.h"

static void
FUNCTION(fft_halfcomplex,pass_2) (const BASE in[],
                                  const size_t istride,
                                  MpIeee out[],
                                  const size_t ostride,
                                  const size_t product,
                                  const size_t n,
                                  const TYPE(gsl_complex) twiddle[]);

static void
FUNCTION(fft_halfcomplex,pass_3) (const BASE in[], 
                                  const size_t istride,
                                  MpIeee out[],
                                  const size_t ostride,
                                  const size_t product,
                                  const size_t n,
                                  const TYPE(gsl_complex) twiddle1[],
                                  const TYPE(gsl_complex) twiddle2[]);

static void
FUNCTION(fft_halfcomplex,pass_4) (const BASE in[],
                                  const size_t istride,
                                  MpIeee out[],
                                  const size_t ostride,
                                  const size_t product,
                                  const size_t n,
                                  const TYPE(gsl_complex) twiddle1[],
                                  const TYPE(gsl_complex) twiddle2[],
                                  const TYPE(gsl_complex) twiddle3[]);

static void
FUNCTION(fft_halfcomplex,pass_5) (const BASE in[],
                                  const size_t istride,
                                  MpIeee out[],
                                  const size_t ostride,
                                  const size_t product,
                                  const size_t n,
                                  const TYPE(gsl_complex) twiddle1[],
                                  const TYPE(gsl_complex) twiddle2[],
                                  const TYPE(gsl_complex) twiddle3[],
                                  const TYPE(gsl_complex) twiddle4[]);

static void
FUNCTION(fft_halfcomplex,pass_n) (const BASE in[],
                                  const size_t istride,
                                  MpIeee out[],
                                  const size_t ostride,
                                  const size_t factor,
                                  const size_t product,
                                  const size_t n,
                                  const TYPE(gsl_complex) twiddle[]);




