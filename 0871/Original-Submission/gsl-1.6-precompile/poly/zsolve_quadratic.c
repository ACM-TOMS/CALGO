#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* poly/zsolve_quadratic.c
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

/* complex_solve_quadratic.c - finds complex roots of a x^2 + b x + c = 0 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_poly.h>

int
 gsl_poly_complex_solve_quadratic(MpIeee a, MpIeee b, MpIeee c,
                                  gsl_complex *z0, gsl_complex *z1)
{
  MpIeee disc=  b * b - MpIeee( "4" ) * a * c;

  if (disc > MpIeee( "0" ))
    {
      if (b == MpIeee( "0" ))
        {
          MpIeee s=  fabs (MpIeee( "0.5" ) * sqrt (disc) / a);
          GSL_REAL (*z0) = -s;
          GSL_IMAG (*z0) = 0;
          GSL_REAL (*z1) = s;
          GSL_IMAG (*z1) = 0;
        }
      else
        {
          MpIeee sgnb=  (b > MpIeee( "0" ) ? MpIeee( "1" ) : -MpIeee( "1" ));
          MpIeee temp=  -MpIeee( "0.5" ) * (b + sgnb * sqrt (disc));
          MpIeee r1=  temp / a;
          MpIeee r2=  c / temp;

          if (r1 < r2)
            {
              GSL_REAL (*z0) = r1;
              GSL_IMAG (*z0) = 0;
              GSL_REAL (*z1) = r2;
              GSL_IMAG (*z1) = 0;
            }
          else
            {
              GSL_REAL (*z0) = r2;
              GSL_IMAG (*z0) = 0;
              GSL_REAL (*z1) = r1;
              GSL_IMAG (*z1) = 0;
            }
        }
      return 2;
    }
  else if (disc == MpIeee( "0" ))
    {
      GSL_REAL (*z0) = -0.5 * b / a;
      GSL_IMAG (*z0) = 0;

      GSL_REAL (*z1) = -0.5 * b / a;
      GSL_IMAG (*z1) = 0;

      return 2;
    }
  else
    {
      MpIeee s=  fabs (MpIeee( "0.5" ) * sqrt (-disc) / a);

      GSL_REAL (*z0) = -0.5 * b / a;
      GSL_IMAG (*z0) = -s;

      GSL_REAL (*z1) = -0.5 * b / a;
      GSL_IMAG (*z1) = s;

      return 2;
    }
}
