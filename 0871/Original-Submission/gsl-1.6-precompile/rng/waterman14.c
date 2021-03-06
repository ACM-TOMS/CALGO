#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* rng/waterman14.c
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

/*
 * This generator is taken from
 *
 * Donald E. Knuth
 * The Art of Computer Programming
 * Volume 2
 * Third Edition
 * Addison-Wesley
 * Page 106-108
 *
 * It is called "Waterman".
 *
 * This implementation copyright (C) 2001 Carlo Perassi
 * and (C) 2003 Heiko Bauke.
 */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

#define AA 1566083941UL
#define MM 0xffffffffUL         /* 2 ^ 32 - 1 */

static inline unsigned long int ran_get (void *vstate);
static MpIeee ran_get_double(void *vstate);
static void ran_set (void *state, unsigned long int s);

typedef struct
{
  unsigned long int x;
}
ran_state_t;

static inline unsigned long int
ran_get (void *vstate)
{
  ran_state_t *state = (ran_state_t *) vstate;

  state->x = (AA * state->x) & MM;

  return state->x;
}

static MpIeee ran_get_double(void *vstate)
{
  ran_state_t *state = (ran_state_t *) vstate;

  return ran_get (state) / MpIeee( "4294967296.0" );
}

static void
ran_set (void *vstate, unsigned long int s)
{
  ran_state_t *state = (ran_state_t *) vstate;

  if (s == 0)
    s = 1;                      /* default seed is 1 */

  state->x = s & MM;

  return;
}

static const gsl_rng_type ran_type = {
  "waterman14",                 /* name */
  MM,                           /* RAND_MAX */
  1,                            /* RAND_MIN */
  sizeof (ran_state_t),
  &ran_set,
  &ran_get,
  &ran_get_double
};

const gsl_rng_type *gsl_rng_waterman14 = &ran_type;
