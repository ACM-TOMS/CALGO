#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* cheb/init.c
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

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_chebyshev.h>

/*-*-*-*-*-*-*-*-*-*-*-* Allocators *-*-*-*-*-*-*-*-*-*-*-*/

gsl_cheb_series * 
gsl_cheb_alloc(const size_t order)
{
  gsl_cheb_series * cs = (gsl_cheb_series *) malloc(sizeof(gsl_cheb_series));
  
  if(cs == 0) {
    GSL_ERROR_VAL("failed to allocate gsl_cheb_series struct", GSL_ENOMEM, 0);
  }
  
  cs->order    = order;
  cs->order_sp = order;

  cs->c = (MpIeee*) malloc((order+1) * sizeof(MpIeee));

  if(cs->c == 0) {
    GSL_ERROR_VAL("failed to allocate cheb coefficients", GSL_ENOMEM, 0);
  }

  cs->f = (MpIeee*) malloc((order+1) * sizeof(MpIeee));

  if(cs->f == 0) {
    GSL_ERROR_VAL("failed to allocate cheb function space", GSL_ENOMEM, 0);
  }

  return cs;
}


void gsl_cheb_free(gsl_cheb_series * cs)
{
  free(cs->f);
  free(cs->c);
  free(cs);
}

/*-*-*-*-*-*-*-*-*-*-*-* Initializer *-*-*-*-*-*-*-*-*-*-*-*/

int  gsl_cheb_init(gsl_cheb_series * cs, const gsl_function *func,
                  const MpIeee a, const MpIeee b)
{
  size_t k, j;

  if(a >= b) {
    GSL_ERROR_VAL("null function interval [a,b]", GSL_EDOM, 0);
  }
  cs->a = a;
  cs->b = b;
  /* cs->err = 0.0; */

  { 
    MpIeee bma=  MpIeee( "0.5" ) * (cs->b - cs->a);
    MpIeee bpa=  MpIeee( "0.5" ) * (cs->b + cs->a);
    MpIeee fac=  MpIeee( "2.0" )/(cs->order +MpIeee( "1.0" ));

    for(k = 0; k<=cs->order; k++) {
      MpIeee y=  cos(M_PI * (k+MpIeee( "0.5" ))/(cs->order+MpIeee( "1" )));
      cs->f[k] = GSL_FN_EVAL(func, (y*bma + bpa));
    }
    
    for(j = 0; j<=cs->order; j++) {
      MpIeee sum=  MpIeee( "0.0" );
      for(k = 0; k<=cs->order; k++) 
        sum += cs->f[k]*cos(M_PI * j*(k+MpIeee( "0.5" ))/(cs->order+MpIeee( "1" )));
      cs->c[j] = fac * sum;
    }
    
  }
  return GSL_SUCCESS;
}





