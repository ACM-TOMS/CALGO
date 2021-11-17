#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* complex/math.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Jorma Olavi Tähtinen, Brian Gough
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

/* Basic complex arithmetic functions

 * Original version by Jorma Olavi Tähtinen <jotahtin@cc.hut.fi>
 *
 * Modified for GSL by Brian Gough, 3/2000
 */

/* The following references describe the methods used in these
 * functions,
 *
 *   T. E. Hull and Thomas F. Fairgrieve and Ping Tak Peter Tang,
 *   "Implementing Complex Elementary Functions Using Exception
 *   Handling", ACM Transactions on Mathematical Software, Volume 20
 *   (1994), pp 215-244, Corrigenda, p553
 *
 *   Hull et al, "Implementing the complex arcsin and arccosine
 *   functions using exception handling", ACM Transactions on
 *   Mathematical Software, Volume 23 (1997) pp 299-335
 *
 *   Abramowitz and Stegun, Handbook of Mathematical Functions, "Inverse
 *   Circular Functions in Terms of Real and Imaginary Parts", Formulas
 *   4.4.37, 4.4.38, 4.4.39
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

/**********************************************************************
 * Complex numbers
 **********************************************************************/

#ifndef HIDE_INLINE_STATIC
gsl_complex
gsl_complex_rect (MpIeee x, MpIeee y)
{                               /* return z = x + i y */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, x, y);
  return z;
}
#endif

gsl_complex
gsl_complex_polar (MpIeee r, MpIeee theta)
{                               /* return z = r exp(i theta) */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, r * cos (theta), r * sin (theta));
  return z;
}

/**********************************************************************
 * Properties of complex numbers
 **********************************************************************/

MpIeee gsl_complex_arg(gsl_complex z)
{                               /* return arg(z),  -pi < arg(z) <= +pi */
  MpIeee x=  GSL_REAL (z);
  MpIeee y=  GSL_IMAG (z);

  if (x == MpIeee( "0.0" ) && y == MpIeee( "0.0" ))
    {
      return MpIeee( "0" );
    }

  return atan(y,x);//atan2 (y, x);
}

MpIeee gsl_complex_abs(gsl_complex z)
{                               /* return |z| */
  return hypot (GSL_REAL (z), GSL_IMAG (z));
}

MpIeee gsl_complex_abs2(gsl_complex z)
{                               /* return |z|^2 */
  MpIeee x=  GSL_REAL (z);
  MpIeee y=  GSL_IMAG (z);

  return (x * x + y * y);
}

MpIeee gsl_complex_logabs(gsl_complex z)
{                               /* return log|z| */
  MpIeee xabs=  fabs (GSL_REAL (z));
  MpIeee yabs=  fabs (GSL_IMAG (z));
  MpIeee max;MpIeee  u;

  if (xabs >= yabs)
    {
      max = xabs;
      u = yabs / xabs;
    }
  else
    {
      max = yabs;
      u = xabs / yabs;
    }

  /* Handle underflow when u is close to 0 */

  return log (max) + MpIeee( "0.5" ) * log1p (u * u);
}


/***********************************************************************
 * Complex arithmetic operators
 ***********************************************************************/

gsl_complex
gsl_complex_add (gsl_complex a, gsl_complex b)
{                               /* z=a+b */
  MpIeee ar=  GSL_REAL (a);MpIeee  ai=  GSL_IMAG (a);
  MpIeee br=  GSL_REAL (b);MpIeee  bi=  GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar + br, ai + bi);
  return z;
}

gsl_complex
gsl_complex_add_real (gsl_complex a, MpIeee x)
{                               /* z=a+x */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a) + x, GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_add_imag (gsl_complex a, MpIeee y)
{                               /* z=a+iy */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a), GSL_IMAG (a) + y);
  return z;
}


gsl_complex
gsl_complex_sub (gsl_complex a, gsl_complex b)
{                               /* z=a-b */
  MpIeee ar=  GSL_REAL (a);MpIeee  ai=  GSL_IMAG (a);
  MpIeee br=  GSL_REAL (b);MpIeee  bi=  GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar - br, ai - bi);
  return z;
}

gsl_complex
gsl_complex_sub_real (gsl_complex a, MpIeee x)
{                               /* z=a-x */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a) - x, GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_sub_imag (gsl_complex a, MpIeee y)
{                               /* z=a-iy */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a), GSL_IMAG (a) - y);
  return z;
}

gsl_complex
gsl_complex_mul (gsl_complex a, gsl_complex b)
{                               /* z=a*b */
  MpIeee ar=  GSL_REAL (a);MpIeee  ai=  GSL_IMAG (a);
  MpIeee br=  GSL_REAL (b);MpIeee  bi=  GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar * br - ai * bi, ar * bi + ai * br);
  return z;
}

gsl_complex
gsl_complex_mul_real (gsl_complex a, MpIeee x)
{                               /* z=a*x */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, x * GSL_REAL (a), x * GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_mul_imag (gsl_complex a, MpIeee y)
{                               /* z=a*iy */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, -y * GSL_IMAG (a), y * GSL_REAL (a));
  return z;
}

gsl_complex
gsl_complex_div (gsl_complex a, gsl_complex b)
{                               /* z=a/b */
  MpIeee ar=  GSL_REAL (a);MpIeee  ai=  GSL_IMAG (a);
  MpIeee br=  GSL_REAL (b);MpIeee  bi=  GSL_IMAG (b);

  MpIeee s=  MpIeee( "1.0" ) / gsl_complex_abs (b);

  MpIeee sbr=  s * br;
  MpIeee sbi=  s * bi;

  MpIeee zr=  (ar * sbr + ai * sbi) * s;
  MpIeee zi=  (ai * sbr - ar * sbi) * s;

  gsl_complex z;
  GSL_SET_COMPLEX (&z, zr, zi);
  return z;
}

gsl_complex
gsl_complex_div_real (gsl_complex a, MpIeee x)
{                               /* z=a/x */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a) / x, GSL_IMAG (a) / x);
  return z;
}

gsl_complex
gsl_complex_div_imag (gsl_complex a, MpIeee y)
{                               /* z=a/(iy) */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_IMAG (a) / y,  - GSL_REAL (a) / y);
  return z;
}

gsl_complex
gsl_complex_conjugate (gsl_complex a)
{                               /* z=conj(a) */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a), -GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_negative (gsl_complex a)
{                               /* z=-a */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, -GSL_REAL (a), -GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_inverse (gsl_complex a)
{                               /* z=1/a */
  MpIeee s=  MpIeee( "1.0" ) / gsl_complex_abs (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, (GSL_REAL (a) * s) * s, -(GSL_IMAG (a) * s) * s);
  return z;
}

/**********************************************************************
 * Elementary complex functions
 **********************************************************************/

gsl_complex
gsl_complex_sqrt (gsl_complex a)
{                               /* z=sqrt(a) */
  gsl_complex z;

  if (GSL_REAL (a) == 0.0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, 0, 0);
    }
  else
    {
      MpIeee x=  fabs (GSL_REAL (a));
      MpIeee y=  fabs (GSL_IMAG (a));
      MpIeee w;

      if (x >= y)
        {
          MpIeee t=  y / x;
          w = sqrt (x) * sqrt (MpIeee( "0.5" ) * (MpIeee( "1.0" ) + sqrt (MpIeee( "1.0" ) + t * t)));
        }
      else
        {
          MpIeee t=  x / y;
          w = sqrt (y) * sqrt (MpIeee( "0.5" ) * (t + sqrt (MpIeee( "1.0" ) + t * t)));
        }

      if (GSL_REAL (a) >= 0.0)
        {
          MpIeee ai=  GSL_IMAG (a);
          GSL_SET_COMPLEX (&z, w, ai / (2.0 * w));
        }
      else
        {
          MpIeee ai=  GSL_IMAG (a);
          MpIeee vi=  (ai >= MpIeee( "0" )) ? w : -w;
          GSL_SET_COMPLEX (&z, ai / (2.0 * vi), vi);
        }
    }

  return z;
}

gsl_complex
gsl_complex_sqrt_real (MpIeee x)
{                               /* z=sqrt(x) */
  gsl_complex z;

  if (x >= MpIeee( "0" ))
    {
      GSL_SET_COMPLEX (&z, sqrt (x), 0.0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, 0.0, sqrt (-x));
    }

  return z;
}

gsl_complex
gsl_complex_exp (gsl_complex a)
{                               /* z=exp(a) */
  MpIeee rho=  exp (GSL_REAL (a));
  MpIeee theta=  GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, rho * cos (theta), rho * sin (theta));
  return z;
}

gsl_complex
gsl_complex_pow (gsl_complex a, gsl_complex b)
{                               /* z=a^b */
  gsl_complex z;

  if (GSL_REAL (a) == 0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, 0.0, 0.0);
    }
  else
    {
      MpIeee logr=  gsl_complex_logabs (a);
      MpIeee theta=  gsl_complex_arg (a);

      MpIeee br=  GSL_REAL (b);MpIeee  bi=  GSL_IMAG (b);

      MpIeee rho=  exp (logr * br - bi * theta);
      MpIeee beta=  theta * br + bi * logr;

      GSL_SET_COMPLEX (&z, rho * cos (beta), rho * sin (beta));
    }

  return z;
}

gsl_complex
gsl_complex_pow_real (gsl_complex a, MpIeee b)
{                               /* z=a^b */
  gsl_complex z;

  if (GSL_REAL (a) == 0 && GSL_IMAG (a) == 0)
    {
      GSL_SET_COMPLEX (&z, 0, 0);
    }
  else
    {
      MpIeee logr=  gsl_complex_logabs (a);
      MpIeee theta=  gsl_complex_arg (a);
      MpIeee rho=  exp (logr * b);
      MpIeee beta=  theta * b;
      GSL_SET_COMPLEX (&z, rho * cos (beta), rho * sin (beta));
    }

  return z;
}

gsl_complex
gsl_complex_log (gsl_complex a)
{                               /* z=log(a) */
  MpIeee logr=  gsl_complex_logabs (a);
  MpIeee theta=  gsl_complex_arg (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, logr, theta);
  return z;
}

gsl_complex
gsl_complex_log10 (gsl_complex a)
{                               /* z = log10(a) */
  return gsl_complex_mul_real (gsl_complex_log (a), 1 / log (10.));
}

gsl_complex
gsl_complex_log_b (gsl_complex a, gsl_complex b)
{
  return gsl_complex_div (gsl_complex_log (a), gsl_complex_log (b));
}

/***********************************************************************
 * Complex trigonometric functions
 ***********************************************************************/

gsl_complex
gsl_complex_sin (gsl_complex a)
{                               /* z = sin(a) */
  MpIeee R=  GSL_REAL (a);MpIeee  I=  GSL_IMAG (a);

  gsl_complex z;

  if (I == MpIeee( "0.0" )) 
    {
      /* avoid returing negative zero (-0.0) for the imaginary part  */

      GSL_SET_COMPLEX (&z, sin (R), 0.0);  
    } 
  else 
    {
      GSL_SET_COMPLEX (&z, sin (R) * cosh (I), cos (R) * sinh (I));
    }

  return z;
}

gsl_complex
gsl_complex_cos (gsl_complex a)
{                               /* z = cos(a) */
  MpIeee R=  GSL_REAL (a);MpIeee  I=  GSL_IMAG (a);

  gsl_complex z;

  if (I == MpIeee( "0.0" )) 
    {
      /* avoid returing negative zero (-0.0) for the imaginary part  */

      GSL_SET_COMPLEX (&z, cos (R), 0.0);  
    } 
  else 
    {
      GSL_SET_COMPLEX (&z, cos (R) * cosh (I), sin (R) * sinh (-I));
    }

  return z;
}

gsl_complex
gsl_complex_tan (gsl_complex a)
{                               /* z = tan(a) */
  MpIeee R=  GSL_REAL (a);MpIeee  I=  GSL_IMAG (a);

  gsl_complex z;

  if (fabs (I) < 1)
    {
      MpIeee D=  pow (cos (R), MpIeee( "2.0" )) + pow (sinh (I), MpIeee( "2.0" ));

      GSL_SET_COMPLEX (&z, 0.5 * sin (2 * R) / D, 0.5 * sinh (2 * I) / D);
    }
  else
    {
      MpIeee u=  exp (-I);
      MpIeee C=  MpIeee( "2" ) * u / (MpIeee( "1" ) - pow (u, MpIeee( "2.0" )));
      MpIeee D=  MpIeee( "1" ) + pow (cos (R), MpIeee( "2.0" )) * pow (C, MpIeee( "2.0" ));

      MpIeee S=  pow (C, MpIeee( "2.0" ));
      MpIeee T=  MpIeee( "1.0" ) / tanh (I);

      GSL_SET_COMPLEX (&z, 0.5 * sin (2 * R) * S / D, T / D);
    }

  return z;
}

gsl_complex
gsl_complex_sec (gsl_complex a)
{                               /* z = sec(a) */
  gsl_complex z = gsl_complex_cos (a);
  return gsl_complex_inverse (z);
}

gsl_complex
gsl_complex_csc (gsl_complex a)
{                               /* z = csc(a) */
  gsl_complex z = gsl_complex_sin (a);
  return gsl_complex_inverse(z);
}


gsl_complex
gsl_complex_cot (gsl_complex a)
{                               /* z = cot(a) */
  gsl_complex z = gsl_complex_tan (a);
  return gsl_complex_inverse (z);
}

/**********************************************************************
 * Inverse Complex Trigonometric Functions
 **********************************************************************/

gsl_complex
gsl_complex_arcsin (gsl_complex a)
{                               /* z = arcsin(a) */
  MpIeee R=  GSL_REAL (a);MpIeee  I=  GSL_IMAG (a);
  gsl_complex z;

  if (I == MpIeee( "0" ))
    {
      z = gsl_complex_arcsin_real (R);
    }
  else
    {
      MpIeee x=  fabs (R);MpIeee  y=  fabs (I);
      MpIeee r=  hypot (x + MpIeee( "1" ), y);MpIeee  s=  hypot (x - MpIeee( "1" ), y);
      MpIeee A=  MpIeee( "0.5" ) * (r + s);
      MpIeee B=  x / A;
      MpIeee y2=  y * y;

      MpIeee real;MpIeee  imag;

      const MpIeee A_crossover=  1.5;const MpIeee  B_crossover=  0.6417;

      if (B <= B_crossover)
        {
          real = asin (B);
        }
      else
        {
          if (x <= MpIeee( "1" ))
            {
              MpIeee D=  MpIeee( "0.5" ) * (A + x) * (y2 / (r + x + MpIeee( "1" )) + (s + (MpIeee( "1" ) - x)));
              real = atan (x / sqrt (D));
            }
          else
            {
              MpIeee Apx=  A + x;
              MpIeee D=  MpIeee( "0.5" ) * (Apx / (r + x + MpIeee( "1" )) + Apx / (s + (x - MpIeee( "1" ))));
              real = atan (x / (y * sqrt (D)));
            }
        }

      if (A <= A_crossover)
        {
          MpIeee Am1;

          if (x < MpIeee( "1" ))
            {
              Am1 = MpIeee( "0.5" ) * (y2 / (r + (x + MpIeee( "1" ))) + y2 / (s + (MpIeee( "1" ) - x)));
            }
          else
            {
              Am1 = MpIeee( "0.5" ) * (y2 / (r + (x + MpIeee( "1" ))) + (s + (x - MpIeee( "1" ))));
            }

          imag = log1p (Am1 + sqrt (Am1 * (A + MpIeee( "1" ))));
        }
      else
        {
          imag = log (A + sqrt (A * A - MpIeee( "1" )));
        }

      GSL_SET_COMPLEX (&z, (R >= MpIeee( "0" )) ? real : -real, (I >= MpIeee( "0" )) ? imag : -imag);
    }

  return z;
}

gsl_complex
gsl_complex_arcsin_real (MpIeee a)
{                               /* z = arcsin(a) */
  gsl_complex z;

  if (fabs (a) <= 1.0)
    {
      GSL_SET_COMPLEX (&z, asin (a), 0.0);
    }
  else
    {
      if (a < MpIeee( "0.0" ))
        {
          GSL_SET_COMPLEX (&z, -M_PI_2, acosh (-a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, M_PI_2, -acosh (a));
        }
    }

  return z;
}

gsl_complex
gsl_complex_arccos (gsl_complex a)
{                               /* z = arccos(a) */
  MpIeee R=  GSL_REAL (a);MpIeee  I=  GSL_IMAG (a);
  gsl_complex z;

  if (I == MpIeee( "0" ))
    {
      z = gsl_complex_arccos_real (R);
    }
  else
    {
      MpIeee x=  fabs (R);MpIeee  y=  fabs (I);
      MpIeee r=  hypot (x + MpIeee( "1" ), y);MpIeee  s=  hypot (x - MpIeee( "1" ), y);
      MpIeee A=  MpIeee( "0.5" ) * (r + s);
      MpIeee B=  x / A;
      MpIeee y2=  y * y;

      MpIeee real;MpIeee  imag;

      const MpIeee A_crossover=  1.5;const MpIeee  B_crossover=  0.6417;

      if (B <= B_crossover)
        {
          real = acos (B);
        }
      else
        {
          if (x <= MpIeee( "1" ))
            {
              MpIeee D=  MpIeee( "0.5" ) * (A + x) * (y2 / (r + x + MpIeee( "1" )) + (s + (MpIeee( "1" ) - x)));
              real = atan (sqrt (D) / x);
            }
          else
            {
              MpIeee Apx=  A + x;
              MpIeee D=  MpIeee( "0.5" ) * (Apx / (r + x + MpIeee( "1" )) + Apx / (s + (x - MpIeee( "1" ))));
              real = atan ((y * sqrt (D)) / x);
            }
        }

      if (A <= A_crossover)
        {
          MpIeee Am1;

          if (x < MpIeee( "1" ))
            {
              Am1 = MpIeee( "0.5" ) * (y2 / (r + (x + MpIeee( "1" ))) + y2 / (s + (MpIeee( "1" ) - x)));
            }
          else
            {
              Am1 = MpIeee( "0.5" ) * (y2 / (r + (x + MpIeee( "1" ))) + (s + (x - MpIeee( "1" ))));
            }

          imag = log1p (Am1 + sqrt (Am1 * (A + MpIeee( "1" ))));
        }
      else
        {
          imag = log (A + sqrt (A * A - MpIeee( "1" )));
        }

      GSL_SET_COMPLEX (&z, (R >= MpIeee( "0" )) ? real : M_PI - real, (I >= MpIeee( "0" )) ? -imag : imag);
    }

  return z;
}

gsl_complex
gsl_complex_arccos_real (MpIeee a)
{                               /* z = arccos(a) */
  gsl_complex z;

  if (fabs (a) <= 1.0)
    {
      GSL_SET_COMPLEX (&z, acos (a), 0);
    }
  else
    {
      if (a < MpIeee( "0.0" ))
        {
          GSL_SET_COMPLEX (&z, M_PI, -acosh (-a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, 0, acosh (a));
        }
    }

  return z;
}

gsl_complex
gsl_complex_arctan (gsl_complex a)
{                               /* z = arctan(a) */
  MpIeee R=  GSL_REAL (a);MpIeee  I=  GSL_IMAG (a);
  gsl_complex z;

  if (I == MpIeee( "0" ))
    {
      GSL_SET_COMPLEX (&z, atan (R), 0);
    }
  else
    {
      /* FIXME: This is a naive implementation which does not fully
         take into account cancellation errors, overflow, underflow
         etc.  It would benefit from the Hull et al treatment. */

      MpIeee r=  hypot (R, I);

      MpIeee imag;

      MpIeee u=  MpIeee( "2" ) * I / (MpIeee( "1" ) + r * r);

      /* FIXME: the following cross-over should be optimized but 0.1
         seems to work ok */

      if (fabs (u) < 0.1)
        {
          imag = MpIeee( "0.25" ) * (log1p (u) - log1p (-u));
        }
      else
        {
          MpIeee A=  hypot (R, I + MpIeee( "1" ));
          MpIeee B=  hypot (R, I - MpIeee( "1" ));
          imag = MpIeee( "0.5" ) * log (A / B);
        }

      if (R == MpIeee( "0" ))
        {
          if (I > MpIeee( "1" ))
            {
              GSL_SET_COMPLEX (&z, M_PI_2, imag);
            }
          else if (I < -MpIeee( "1" ))
            {
              GSL_SET_COMPLEX (&z, -M_PI_2, imag);
            }
          else
            {
              GSL_SET_COMPLEX (&z, 0, imag);
            };
        }
      else
        {
          //GSL_SET_COMPLEX (&z, 0.5 * atan (2 * R, ((1 + r) * (1 - r)))/ imag); //changed atan2 (x,y) into atan(x/y) for MpIeee...
        }
    }

  return z;
}

gsl_complex
gsl_complex_arcsec (gsl_complex a)
{                               /* z = arcsec(a) */
  gsl_complex z = gsl_complex_inverse (a);
  return gsl_complex_arccos (z);
}

gsl_complex
gsl_complex_arcsec_real (MpIeee a)
{                               /* z = arcsec(a) */
  gsl_complex z;

  if (a <= -MpIeee( "1.0" ) || a >= MpIeee( "1.0" ))
    {
      GSL_SET_COMPLEX (&z, acos (1 / a), 0.0);
    }
  else
    {
      if (a >= MpIeee( "0.0" ))
        {
          GSL_SET_COMPLEX (&z, 0, acosh (1 / a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, M_PI, -acosh (-1 / a));
        }
    }

  return z;
}

gsl_complex
gsl_complex_arccsc (gsl_complex a)
{                               /* z = arccsc(a) */
  gsl_complex z = gsl_complex_inverse (a);
  return gsl_complex_arcsin (z);
}

gsl_complex
gsl_complex_arccsc_real (MpIeee a)
{                               /* z = arccsc(a) */
  gsl_complex z;

  if (a <= -MpIeee( "1.0" ) || a >= MpIeee( "1.0" ))
    {
      GSL_SET_COMPLEX (&z, asin (1 / a), 0.0);
    }
  else
    {
      if (a >= MpIeee( "0.0" ))
        {
          GSL_SET_COMPLEX (&z, M_PI_2, -acosh (1 / a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, -M_PI_2, acosh (-1 / a));
        }
    }

  return z;
}

gsl_complex
gsl_complex_arccot (gsl_complex a)
{                               /* z = arccot(a) */
  gsl_complex z;

  if (GSL_REAL (a) == 0.0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, M_PI_2, 0);
    }
  else
    {
      z = gsl_complex_inverse (a);
      z = gsl_complex_arctan (z);
    }

  return z;
}

/**********************************************************************
 * Complex Hyperbolic Functions
 **********************************************************************/

gsl_complex
gsl_complex_sinh (gsl_complex a)
{                               /* z = sinh(a) */
  MpIeee R=  GSL_REAL (a);MpIeee  I=  GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, sinh (R) * cos (I), cosh (R) * sin (I));
  return z;
}

gsl_complex
gsl_complex_cosh (gsl_complex a)
{                               /* z = cosh(a) */
  MpIeee R=  GSL_REAL (a);MpIeee  I=  GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, cosh (R) * cos (I), sinh (R) * sin (I));
  return z;
}

gsl_complex
gsl_complex_tanh (gsl_complex a)
{                               /* z = tanh(a) */
  MpIeee R=  GSL_REAL (a);MpIeee  I=  GSL_IMAG (a);

  gsl_complex z;

  if (fabs(R) < 1.0) 
    {
      MpIeee D=  pow (cos (I), MpIeee( "2.0" )) + pow (sinh (R), MpIeee( "2.0" ));
      
      GSL_SET_COMPLEX (&z, sinh (R) * cosh (R) / D, 0.5 * sin (2 * I) / D);
    }
  else
    {
      MpIeee D=  pow (cos (I), MpIeee( "2.0" )) + pow (sinh (R), MpIeee( "2.0" ));
      MpIeee F=  MpIeee( "1" ) + pow (cos (I) / sinh (R), MpIeee( "2.0" ));

      GSL_SET_COMPLEX (&z, 1.0 / (tanh (R) * F), 0.5 * sin (2 * I) / D);
    }

  return z;
}

gsl_complex
gsl_complex_sech (gsl_complex a)
{                               /* z = sech(a) */
  gsl_complex z = gsl_complex_cosh (a);
  return gsl_complex_inverse (z);
}

gsl_complex
gsl_complex_csch (gsl_complex a)
{                               /* z = csch(a) */
  gsl_complex z = gsl_complex_sinh (a);
  return gsl_complex_inverse (z);
}

gsl_complex
gsl_complex_coth (gsl_complex a)
{                               /* z = coth(a) */
  gsl_complex z = gsl_complex_tanh (a);
  return gsl_complex_inverse (z);
}

/**********************************************************************
 * Inverse Complex Hyperbolic Functions
 **********************************************************************/

gsl_complex
gsl_complex_arcsinh (gsl_complex a)
{                               /* z = arcsinh(a) */
  gsl_complex z = gsl_complex_mul_imag(a, 1.0);
  z = gsl_complex_arcsin (z);
  z = gsl_complex_mul_imag (z, -1.0);
  return z;
}

gsl_complex
gsl_complex_arccosh (gsl_complex a)
{                               /* z = arccosh(a) */
  gsl_complex z = gsl_complex_arccos (a);
  z = gsl_complex_mul_imag (z, GSL_IMAG(z) > 0 ? -1.0 : 1.0);
  return z;
}

gsl_complex
gsl_complex_arccosh_real (MpIeee a)
{                               /* z = arccosh(a) */
  gsl_complex z;

  if (a >= MpIeee( "1" ))
    {
      GSL_SET_COMPLEX (&z, acosh (a), 0);
    }
  else
    {
      if (a >= -MpIeee( "1.0" ))
        {
          GSL_SET_COMPLEX (&z, 0, acos (a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, acosh (-a), M_PI);
        }
    }

  return z;
}

gsl_complex
gsl_complex_arctanh (gsl_complex a)
{                               /* z = arctanh(a) */
  if (GSL_IMAG (a) == 0.0)
    {
      return gsl_complex_arctanh_real (GSL_REAL (a));
    }
  else
    {
      gsl_complex z = gsl_complex_mul_imag(a, 1.0);
      z = gsl_complex_arctan (z);
      z = gsl_complex_mul_imag (z, -1.0);
      return z;
    }
}

gsl_complex
gsl_complex_arctanh_real (MpIeee a)
{                               /* z = arctanh(a) */
  gsl_complex z;

  if (a > -MpIeee( "1.0" ) && a < MpIeee( "1.0" ))
    {
      GSL_SET_COMPLEX (&z, atanh (a), 0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, atanh (1 / a), (a < MpIeee( "0" )) ? M_PI_2 : -M_PI_2);
    }

  return z;
}

gsl_complex
gsl_complex_arcsech (gsl_complex a)
{                               /* z = arcsech(a); */
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arccosh (t);
}

gsl_complex
gsl_complex_arccsch (gsl_complex a)
{                               /* z = arccsch(a) */
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arcsinh (t);
}

gsl_complex
gsl_complex_arccoth (gsl_complex a)
{                               /* z = arccoth(a) */
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arctanh (t);
}
