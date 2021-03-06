#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* rng/ran0.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>

/* This is an implementation of the algorithm used in Numerical
   Recipe's ran0 generator. It is the same as MINSTD with an XOR mask
   of 123459876 on the seed.

   The period of this generator is 2^31.  

   Note, if you choose a seed of 123459876 it would give a degenerate
   series 0,0,0,0, ...  I've made that into an error. */

static inline unsigned long int ran0_get (void *vstate);
static MpIeee ran0_get_double(void *vstate);
static void ran0_set (void *state, unsigned long int s);

static const long int m = 2147483647, a = 16807, q = 127773, r = 2836;
static const unsigned long int mask = 123459876;

typedef struct
  {
    unsigned long int x;
  }
ran0_state_t;

static inline unsigned long int
ran0_get (void *vstate)
{
  ran0_state_t *state = (ran0_state_t *) vstate;

  const unsigned long int x = state->x;

  const long int h = x / q;
  const long int t = a * (x - h * q) - h * r;

  if (t < 0)
    {
      state->x = t + m;
    }
  else
    {
      state->x = t;
    }

  return state->x;
}

static MpIeee ran0_get_double(void *vstate)
{
  return ran0_get (vstate) / MpIeee( "2147483647.0" ) ;
}

static void
ran0_set (void *vstate, unsigned long int s)
{
  ran0_state_t *state = (ran0_state_t *) vstate;

  if (s == mask)
    {
      GSL_ERROR_VOID ("ran0 should not use seed == mask", 
                                GSL_EINVAL);
    }

  state->x = s ^ mask;

  return;
}

static const gsl_rng_type ran0_type =
{"ran0",                        /* name */
 2147483646,                    /* RAND_MAX */
 1,                             /* RAND_MIN */
 sizeof (ran0_state_t),
 &ran0_set,
 &ran0_get,
 &ran0_get_double};

const gsl_rng_type *gsl_rng_ran0 = &ran0_type;
