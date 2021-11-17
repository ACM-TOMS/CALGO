#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/beta_inc.c
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
/* Modified for cdfs by Brian Gough, June 2003 */

static MpIeee beta_cont_frac(const MpIeee a, const MpIeee b, const MpIeee x,
                const MpIeee epsabs)
{
  const unsigned int  max_iter=  512;    /* control iterations      */
  const MpIeee cutoff=  2.0 * GSL_DBL_MIN;      /* control the zero cutoff */
  unsigned int  iter_count=  0;
  MpIeee cf;

  /* standard initialization for continued fraction */
  MpIeee num_term=  MpIeee( "1.0" );
  MpIeee den_term=  MpIeee( "1.0" ) - (a + b) * x / (a + MpIeee( "1.0" ));

  if (fabs (den_term) < cutoff)
    den_term = GSL_NAN;

  den_term = MpIeee( "1.0" ) / den_term;
  cf = den_term;

  while (iter_count < max_iter)
    {
      const int k = iter_count + 1;
      MpIeee coeff=  k * (b - k) * x / (((a - MpIeee( "1.0" )) + MpIeee( "2" ) * k) * (a + MpIeee( "2" ) * k));
      MpIeee delta_frac;

      /* first step */
      den_term = MpIeee( "1.0" ) + coeff * den_term;
      num_term = MpIeee( "1.0" ) + coeff / num_term;

      if (fabs (den_term) < cutoff)
        den_term = GSL_NAN;

      if (fabs (num_term) < cutoff)
        num_term = GSL_NAN;

      den_term = MpIeee( "1.0" ) / den_term;

      delta_frac = den_term * num_term;
      cf *= delta_frac;

      coeff = -(a + k) * (a + b + k) * x / ((a + MpIeee( "2" ) * k) * (a + MpIeee( "2" ) * k + MpIeee( "1.0" )));

      /* second step */
      den_term = MpIeee( "1.0" ) + coeff * den_term;
      num_term = MpIeee( "1.0" ) + coeff / num_term;

      if (fabs (den_term) < cutoff)
        den_term = GSL_NAN;

      if (fabs (num_term) < cutoff)
        num_term = GSL_NAN;

      den_term = MpIeee( "1.0" ) / den_term;

      delta_frac = den_term * num_term;
      cf *= delta_frac;

      if (fabs (delta_frac - 1.0) < 2.0 * GSL_DBL_EPSILON)
        break;

      if (cf * fabs (delta_frac - 1.0) < epsabs)
        break;

      ++iter_count;
    }

  if (iter_count >= max_iter)
    return GSL_NAN;

  return cf;
}

/* The function beta_inc_AXPY(A,Y,a,b,x) computes A * beta_inc(a,b,x)
   + Y taking account of possible cancellations when using the
   hypergeometric transformation beta_inc(a,b,x)=1-beta_inc(b,a,1-x).

   It also adjusts the accuracy of beta_inc() to fit the overall
   absolute error when A*beta_inc is added to Y. (e.g. if Y >>
   A*beta_inc then the accuracy of beta_inc can be reduced) */

static MpIeee beta_inc_AXPY(const MpIeee A, const MpIeee Y,
               const MpIeee a, const MpIeee b, const MpIeee x)
{
  if (x == 0.0)
    {
      return A * MpIeee( "0" ) + Y;
    }
  else if (x == 1.0)
    {
      return A * MpIeee( "1" ) + Y;
    }
  else
    {
      MpIeee ln_beta=  gsl_sf_lnbeta (a, b);
      MpIeee ln_pre=  -ln_beta + a * log (x) + b * log1p (-x);

      MpIeee prefactor=  exp (ln_pre);

      if (x < (a + 1.0) / (a + b + 2.0))
        {
          /* Apply continued fraction directly. */
          MpIeee epsabs=  fabs (Y / (A * prefactor / a)) * GSL_DBL_EPSILON;

          MpIeee cf=  beta_cont_frac (a, b, x, epsabs);

          return A * (prefactor * cf / a) + Y;
        }
      else
        {
          /* Apply continued fraction after hypergeometric transformation. */
          MpIeee epsabs= 
            fabs ((A + Y) / (A * prefactor / b)) * GSL_DBL_EPSILON;
          MpIeee cf=  beta_cont_frac (b, a, MpIeee( "1.0" ) - x, epsabs);
          MpIeee term=  prefactor * cf / b;

          if (A == -Y)
            {
              return -A * term;
            }
          else
            {
              return A * (MpIeee( "1" ) - term) + Y;
            }
        }
    }
}

/* Direct series evaluation for testing purposes only */

#if 0
static MpIeee beta_series(const MpIeee a, const MpIeee b, const MpIeee x,
             const MpIeee epsabs)
{
  MpIeee f=  x / (MpIeee( "1" ) - x);
  MpIeee c=  (b - MpIeee( "1" )) / (a + MpIeee( "1" )) * f;
  MpIeee s=  MpIeee( "1" );
  MpIeee n=  MpIeee( "0" );

  s += c;

  do
    {
      n++;
      c *= -f * (MpIeee( "2" ) + n - b) / (MpIeee( "2" ) + n + a);
      s += c;
    }
  while (n < MpIeee( "512" ) && fabs (c) > GSL_DBL_EPSILON * fabs (s) + epsabs);

  s /= (MpIeee( "1" ) - x);

  return s;
}
#endif
