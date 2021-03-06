#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* poly/zsolve_init.c
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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_errno.h>

gsl_poly_complex_workspace * 
gsl_poly_complex_workspace_alloc (size_t n)
{
  size_t nc ;

  gsl_poly_complex_workspace * w ;
  
  if (n == 0)
    {
      GSL_ERROR_VAL ("matrix size n must be positive integer", GSL_EDOM, 0);
    }

  w = (gsl_poly_complex_workspace *) 
    malloc (sizeof(gsl_poly_complex_workspace));

  if (w == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for struct", GSL_ENOMEM, 0);
    }

  nc = n - 1;

  w->nc = nc;

  w->matrix = (MpIeee*) malloc (nc * nc * sizeof(MpIeee));

  if (w->matrix == 0)
    {
      free (w) ;       /* error in constructor, avoid memory leak */
      
      GSL_ERROR_VAL ("failed to allocate space for workspace matrix", 
                        GSL_ENOMEM, 0);
    }

  return w ;
}

void 
gsl_poly_complex_workspace_free (gsl_poly_complex_workspace * w)
{
  free(w->matrix) ;
  free(w);
}
