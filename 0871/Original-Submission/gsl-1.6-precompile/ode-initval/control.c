#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* ode-initval/control.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

gsl_odeiv_control *
gsl_odeiv_control_alloc(const gsl_odeiv_control_type * T)
{
  gsl_odeiv_control * c = 
    (gsl_odeiv_control *) malloc(sizeof(gsl_odeiv_control));

  if(c == 0) 
    {
      GSL_ERROR_NULL ("failed to allocate space for control struct", 
                      GSL_ENOMEM);
    };

  c->type = T;
  c->state = c->type->alloc();

  if (c->state == 0)
    {
      free (c);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_NULL ("failed to allocate space for control state", 
                      GSL_ENOMEM);
    };

  return c;
}

int
 gsl_odeiv_control_init(gsl_odeiv_control * c, 
                       MpIeee eps_abs, MpIeee eps_rel, 
                       MpIeee a_y, MpIeee a_dydt)
{
  return c->type->init (c->state, eps_abs, eps_rel, a_y, a_dydt);
}

void
gsl_odeiv_control_free(gsl_odeiv_control * c)
{
  c->type->free(c->state);
  free(c);
}

const char *
gsl_odeiv_control_name(const gsl_odeiv_control * c)
{
  return c->type->name;
}

int
 gsl_odeiv_control_hadjust(gsl_odeiv_control * c, gsl_odeiv_step * s, const MpIeee y0[], const MpIeee yerr[], const MpIeee dydt[], MpIeee * h)
{
  return c->type->hadjust(c->state, s->dimension, s->type->order(s->state),
                          y0, yerr, dydt, h);
}
