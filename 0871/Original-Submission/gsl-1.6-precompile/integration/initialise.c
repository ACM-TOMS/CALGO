#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* integration/initialise.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001 Brian Gough
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

static inline
void initialise (gsl_integration_workspace * workspace, MpIeee a, MpIeee b);

static inline
void initialise (gsl_integration_workspace * workspace, MpIeee a, MpIeee b)
{
  workspace->size = 0;
  workspace->nrmax = 0;
  workspace->i = 0;
  workspace->alist[0] = a;
  workspace->blist[0] = b;
  workspace->rlist[0] = 0.0;
  workspace->elist[0] = 0.0;
  workspace->order[0] = 0;
  workspace->level[0] = 0;

  workspace->maximum_level = 0;
}
