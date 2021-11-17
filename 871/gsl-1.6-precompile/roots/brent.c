#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* roots/brent.c
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

/* brent.c -- brent root finding algorithm */

#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "roots.h"


typedef struct
  {
    MpIeee a;MpIeee  b;MpIeee  c;MpIeee  d;MpIeee  e;
    MpIeee fa;MpIeee  fb;MpIeee  fc;
  }
brent_state_t;

static int  brent_init(void * vstate, gsl_function * f, MpIeee * root, MpIeee x_lower, MpIeee x_upper);
static int  brent_iterate(void * vstate, gsl_function * f, MpIeee * root, MpIeee * x_lower, MpIeee * x_upper);


static int
 brent_init(void * vstate, gsl_function * f, MpIeee * root, MpIeee x_lower, MpIeee x_upper)
{
  brent_state_t * state = (brent_state_t *) vstate;

  MpIeee f_lower;MpIeee  f_upper;

  *root = MpIeee( "0.5" ) * (x_lower + x_upper) ;

  SAFE_FUNC_CALL (f, x_lower, &f_lower);
  SAFE_FUNC_CALL (f, x_upper, &f_upper);
  
  state->a = x_lower;
  state->fa = f_lower;

  state->b = x_upper;
  state->fb = f_upper;

  state->c = x_upper;
  state->fc = f_upper;

  state->d = x_upper - x_lower ;
  state->e = x_upper - x_lower ;

  if ((f_lower < MpIeee( "0.0" ) && f_upper < MpIeee( "0.0" )) || (f_lower > MpIeee( "0.0" ) && f_upper > MpIeee( "0.0" )))
    {
      GSL_ERROR ("endpoints do not straddle y=0", GSL_EINVAL);
    }

  return GSL_SUCCESS;

}

static int
 brent_iterate(void * vstate, gsl_function * f, MpIeee * root, MpIeee * x_lower, MpIeee * x_upper)
{
  brent_state_t * state = (brent_state_t *) vstate;

  MpIeee tol;MpIeee  m;

  int  ac_equal=  0;

  MpIeee a=  state->a;MpIeee  b=  state->b;MpIeee  c=  state->c;
  MpIeee fa=  state->fa;MpIeee  fb=  state->fb;MpIeee  fc=  state->fc;
  MpIeee d=  state->d;MpIeee  e=  state->e;
  
  if ((fb < MpIeee( "0" ) && fc < MpIeee( "0" )) || (fb > MpIeee( "0" ) && fc > MpIeee( "0" )))
    {
      ac_equal = 1;
      c = a;
      fc = fa;
      d = b - a;
      e = b - a;
    }
  
  if (fabs (fc) < fabs (fb))
    {
      ac_equal = 1;
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
  
  tol = MpIeee( "0.5" ) * MpIeee( GSL_DBL_EPSILON ) * b.fabs();
  m = MpIeee( "0.5" ) * (c - b);
  
  if (fb == MpIeee( "0" ))
    {
      *root = b;
      *x_lower = b;
      *x_upper = b;
      
      return GSL_SUCCESS;
    }
  
  if (fabs (m) <= tol)
    {
      *root = b;

      if (b < c) 
        {
          *x_lower = b;
          *x_upper = c;
        }
      else
        {
          *x_lower = c;
          *x_upper = b;
        }

      return GSL_SUCCESS;
    }
  
  if (fabs (e) < tol || fabs (fa) <= fabs (fb))
    {
      d = m;            /* use bisection */
      e = m;
    }
  else
    {
      MpIeee p;MpIeee  q;MpIeee  r;   /* use inverse cubic interpolation */
      MpIeee s=  fb / fa;
      
      if (ac_equal)
        {
          p = MpIeee( "2" ) * m * s;
          q = MpIeee( "1" ) - s;
        }
      else
        {
          q = fa / fc;
          r = fb / fc;
          p = s * (MpIeee( "2" ) * m * q * (q - r) - (b - a) * (r - MpIeee( "1" )));
          q = (q - MpIeee( "1" )) * (r - MpIeee( "1" )) * (s - MpIeee( "1" ));
        }
      
      if (p > MpIeee( "0" ))
        {
          q = -q;
        }
      else
        {
          p = -p;
        }
      
      if (2 * p < GSL_MIN (MpIeee( "3" ) * m * q - fabs (tol * q), fabs (e * q)))
        {
          e = d;
          d = p / q;
        }
      else
        {
          /* interpolation failed, fall back to bisection */
          
          d = m;
          e = m;
        }
    }
  
  a = b;
  fa = fb;
  
  if (fabs (d) > tol)
    {
      b += d;
    }
  else
    {
      b += (m > MpIeee( "0" ) ? tol : -tol);
    }
  
  SAFE_FUNC_CALL (f, b, &fb);

  state->a = a ;
  state->b = b ;
  state->c = c ;
  state->d = d ;
  state->e = e ;
  state->fa = fa ;
  state->fb = fb ;
  state->fc = fc ;
  
  /* Update the best estimate of the root and bounds on each
     iteration */
  
  *root = b;
  
  if ((fb < MpIeee( "0" ) && fc < MpIeee( "0" )) || (fb > MpIeee( "0" ) && fc > MpIeee( "0" ))) 
    {
      c = a;
    }

  if (b < c)
    {
      *x_lower = b;
      *x_upper = c;
    }
  else
    {
      *x_lower = c;
      *x_upper = b;
    }

  return GSL_SUCCESS ;
}

  
static const gsl_root_fsolver_type brent_type =
{"brent",                               /* name */
 sizeof (brent_state_t),
 &brent_init,
 &brent_iterate};

const gsl_root_fsolver_type  * gsl_root_fsolver_brent = &brent_type;
