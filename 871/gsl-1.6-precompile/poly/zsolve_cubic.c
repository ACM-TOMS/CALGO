#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* poly/zsolve_cubic.c
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

/* zsolve_cubic.c - finds the complex roots of x^3 + a x^2 + b x + c = 0 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_poly.h>

#define SWAP(a,b) do { MpIeee tmp = b ; b = a ; a = tmp ; } while(0)

int
 gsl_poly_complex_solve_cubic(MpIeee a, MpIeee b, MpIeee c, 
                              gsl_complex *z0, gsl_complex *z1, 
                              gsl_complex *z2)
{
  MpIeee q=  (a * a - MpIeee( "3" ) * b);
  MpIeee r=  (MpIeee( "2" ) * a * a * a - MpIeee( "9" ) * a * b + MpIeee( "27" ) * c);

  MpIeee Q=  q / MpIeee( "9" );
  MpIeee R=  r / MpIeee( "54" );

  MpIeee Q3=  Q * Q * Q;
  MpIeee R2=  R * R;

  MpIeee CR2=  MpIeee( "729" ) * r * r;
  MpIeee CQ3=  MpIeee( "2916" ) * q * q * q;

  if (R == MpIeee( "0" ) && Q == MpIeee( "0" ))
    {
      GSL_REAL (*z0) = -a / 3;
      GSL_IMAG (*z0) = 0;
      GSL_REAL (*z1) = -a / 3;
      GSL_IMAG (*z1) = 0;
      GSL_REAL (*z2) = -a / 3;
      GSL_IMAG (*z2) = 0;
      return 3;
    }
  else if (CR2 == CQ3) 
    {
      /* this test is actually R2 == Q3, written in a form suitable
         for exact computation with integers */

      /* Due to finite precision some double roots may be missed, and
         will be considered to be a pair of complex roots z = x +/-
         epsilon i close to the real axis. */

      MpIeee sqrtQ=  sqrt (Q);

      if (R > MpIeee( "0" ))
        {
          GSL_REAL (*z0) = -2 * sqrtQ - a / 3;
          GSL_IMAG (*z0) = 0;
          GSL_REAL (*z1) = sqrtQ - a / 3;
          GSL_IMAG (*z1) = 0;
          GSL_REAL (*z2) = sqrtQ - a / 3;
          GSL_IMAG (*z2) = 0;
        }
      else
        {
          GSL_REAL (*z0) = -sqrtQ - a / 3;
          GSL_IMAG (*z0) = 0;
          GSL_REAL (*z1) = -sqrtQ - a / 3;
          GSL_IMAG (*z1) = 0;
          GSL_REAL (*z2) = 2 * sqrtQ - a / 3;
          GSL_IMAG (*z2) = 0;
        }
      return 3;
    }
  else if (CR2 < CQ3)  /* equivalent to R2 < Q3 */
    {
      MpIeee sqrtQ=  sqrt (Q);
      MpIeee sqrtQ3=  sqrtQ * sqrtQ * sqrtQ;
      MpIeee theta=  acos (R / sqrtQ3);
      MpIeee norm=  -MpIeee( "2" ) * sqrtQ;
      MpIeee r0=  norm * cos (theta / MpIeee( "3" )) - a / MpIeee( "3" );
      MpIeee r1=  norm * cos ((theta + MpIeee( "2.0" ) * M_PI) / MpIeee( "3" )) - a / MpIeee( "3" );
      MpIeee r2=  norm * cos ((theta - MpIeee( "2.0" ) * M_PI) / MpIeee( "3" )) - a / MpIeee( "3" );

      /* Sort r0, r1, r2 into increasing order */

      if (r0 > r1)
        SWAP (r0, r1);

      if (r1 > r2)
        {
          SWAP (r1, r2);

          if (r0 > r1)
            SWAP (r0, r1);
        }

      GSL_REAL (*z0) = r0;
      GSL_IMAG (*z0) = 0;

      GSL_REAL (*z1) = r1;
      GSL_IMAG (*z1) = 0;

      GSL_REAL (*z2) = r2;
      GSL_IMAG (*z2) = 0;

      return 3;
    }
  else
    {
      MpIeee sgnR=  (R >= MpIeee( "0" ) ? MpIeee( "1" ) : -MpIeee( "1" ));
      MpIeee A=  -sgnR * pow (fabs (R) + sqrt (R2 - Q3), MpIeee( "1.0" ) / MpIeee( "3.0" ));
      MpIeee B=  Q / A;

      if (A + B < MpIeee( "0" ))
        {
          GSL_REAL (*z0) = A + B - a / 3;
          GSL_IMAG (*z0) = 0;

          GSL_REAL (*z1) = -0.5 * (A + B) - a / 3;
          GSL_IMAG (*z1) = -(sqrt (3.0) / 2.0) * fabs(A - B);

          GSL_REAL (*z2) = -0.5 * (A + B) - a / 3;
          GSL_IMAG (*z2) = (sqrt (3.0) / 2.0) * fabs(A - B);
        }
      else
        {
          GSL_REAL (*z0) = -0.5 * (A + B) - a / 3;
          GSL_IMAG (*z0) = -(sqrt (3.0) / 2.0) * fabs(A - B);

          GSL_REAL (*z1) = -0.5 * (A + B) - a / 3;
          GSL_IMAG (*z1) = (sqrt (3.0) / 2.0) * fabs(A - B);

          GSL_REAL (*z2) = A + B - a / 3;
          GSL_IMAG (*z2) = 0;
        }

      return 3;
    }
}
