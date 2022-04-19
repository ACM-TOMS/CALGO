#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/bigauss.c
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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* The Bivariate Gaussian probability distribution is 

   p(x,y) dxdy = (1/(2 pi sigma_x sigma_y sqrt(r))) 
                    exp(-(x^2 + y^2 - 2 r x y)/(2c)) dxdy     

*/

void
gsl_ran_bivariate_gaussian (const gsl_rng * r, 
                            MpIeee sigma_x, MpIeee sigma_y, MpIeee rho,
                            MpIeee *x, MpIeee *y)
{
  MpIeee u;MpIeee  v;MpIeee  r2;MpIeee  scale;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      u = -MpIeee( "1" ) + MpIeee( "2" ) * gsl_rng_uniform (r);
      v = -MpIeee( "1" ) + MpIeee( "2" ) * gsl_rng_uniform (r);

      /* see if it is in the unit circle */
      r2 = u * u + v * v;
    }
  while (r2 > MpIeee( "1.0" ) || r2 == MpIeee( "0" ));

  scale = sqrt (-MpIeee( "2.0" ) * log (r2) / r2);

  *x = sigma_x * u * scale;
  *y = sigma_y * (rho * u + sqrt(MpIeee( "1" ) - rho*rho) * v) * scale;
}

MpIeee gsl_ran_bivariate_gaussian_pdf(const MpIeee x, const MpIeee y, 
                                const MpIeee sigma_x, const MpIeee sigma_y,
                                const MpIeee rho)
{
  MpIeee u=  x / sigma_x ;
  MpIeee v=  y / sigma_y ;
  MpIeee c=  MpIeee( "1" ) - rho*rho ;
  MpIeee p=  (MpIeee( "1" ) / (MpIeee( "2" ) * M_PI * sigma_x * sigma_y * sqrt(c))) 
    * exp (-(u * u - MpIeee( "2" ) * rho * u * v + v * v) / (MpIeee( "2" ) * c));
  return p;
}
