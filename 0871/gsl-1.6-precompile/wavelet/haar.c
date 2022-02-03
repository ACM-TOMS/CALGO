#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* wavelet/haar.c
 * 
 * Copyright (C) 2004 Ivo Alxneit
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_wavelet.h>

static const MpIeee ch_2[2] =  { M_SQRT1_2, M_SQRT1_2 };
static const MpIeee cg_2[2] =  { M_SQRT1_2, -(M_SQRT1_2) };

static int
 haar_init(const MpIeee **h1, const MpIeee **g1, const MpIeee **h2,
           const MpIeee **g2, size_t * nc, size_t * offset,
           const size_t member)
{
  if (member != 2)
    {
      return GSL_FAILURE;
    }

  *h1 = ch_2;
  *g1 = cg_2;
  *h2 = ch_2;
  *g2 = cg_2;

  *nc = 2;
  *offset = 0;

  return GSL_SUCCESS;
}

static int
 haar_centered_init(const MpIeee **h1, const MpIeee **g1, const MpIeee **h2,
                    const MpIeee **g2, size_t * nc, size_t * offset,
                    const size_t member)
{
  if (member != 2)
    {
      return GSL_FAILURE;
    }

  *h1 = ch_2;
  *g1 = cg_2;
  *h2 = ch_2;
  *g2 = cg_2;

  *nc = 2;
  *offset = 1;

  return GSL_SUCCESS;
}

static const gsl_wavelet_type haar_type = {
  "haar",
  &haar_init
};

static const gsl_wavelet_type haar_centered_type = {
  "haar-centered",
  &haar_centered_init
};

const gsl_wavelet_type *gsl_wavelet_haar = &haar_type;
const gsl_wavelet_type *gsl_wavelet_haar_centered = &haar_centered_type;
