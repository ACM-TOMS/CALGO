#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* cdf/gauss.c
 *
 * Copyright (C) 2002, 2004 Jason H. Stover.
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA.
 */

/*
 * Computes the cumulative distribution function for the Gaussian
 * distribution using a rational function approximation.  The
 * computation is for the standard Normal distribution, i.e., mean 0
 * and standard deviation 1. If you want to compute Pr(X < t) for a
 * Gaussian random variable X with non-zero mean m and standard
 * deviation sd not equal to 1, find gsl_cdf_ugaussian ((t-m)/sd).
 * This approximation is accurate to at least double precision. The
 * accuracy was verified with a pari-gp script.  The largest error
 * found was about 1.4E-20. The coefficients were derived by Cody.
 *
 * References:
 *
 * W.J. Cody. "Rational Chebyshev Approximations for the Error
 * Function," Mathematics of Computation, v23 n107 1969, 631-637.
 *
 * W. Fraser, J.F Hart. "On the Computation of Rational Approximations
 * to Continuous Functions," Communications of the ACM, v5 1962.
 *
 * W.J. Kennedy Jr., J.E. Gentle. "Statistical Computing." Marcel Dekker. 1980.
 * 
 *  
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>

#ifndef M_1_SQRT2PI
#define M_1_SQRT2PI (M_2_SQRTPI * M_SQRT1_2 / 2.0)
#endif

#define SQRT32 (4.0 * M_SQRT2)

/*
 * IEEE double precision dependent constants.
 *
 * GAUSS_EPSILON: Smallest positive value such that 
 *                gsl_cdf_gaussian(x) > 0.5.
 * GAUSS_XUPPER: Largest value x such that gsl_cdf_gaussian(x) < 1.0.
 * GAUSS_XLOWER: Smallest value x such that gsl_cdf_gaussian(x) > 0.0.
 */

#define GAUSS_EPSILON  (GSL_DBL_EPSILON / 2)
#define GAUSS_XUPPER (8.572)
#define GAUSS_XLOWER (-37.519)

#define GAUSS_SCALE (16.0)

static MpIeee get_del(MpIeee x, MpIeee rational)
{
  MpIeee xsq=  MpIeee( "0.0" );
  MpIeee del=  MpIeee( "0.0" );
  MpIeee result=  MpIeee( "0.0" );

  xsq = floor (x * GAUSS_SCALE) / GAUSS_SCALE;
  del = (x - xsq) * (x + xsq);
  del *= MpIeee( "0.5" );

  result = exp (-MpIeee( "0.5" ) * xsq * xsq) * exp (-MpIeee( "1.0" ) * del) * rational;

  return result;
}

/*
 * Normal cdf for fabs(x) < 0.66291
 */
static MpIeee gauss_small(const MpIeee x)
{
  unsigned int  i;
  MpIeee result=  MpIeee( "0.0" );
  MpIeee xsq;
  MpIeee xnum;
  MpIeee xden;

  const MpIeee a[5] =  {
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113
  };
  const MpIeee b[4] =  {
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956
  };

  xsq = x * x;
  xnum = a[4] * xsq;
  xden = xsq;

  for (i = 0; i < 3; i++)
    {
      xnum = (xnum + a[i]) * xsq;
      xden = (xden + b[i]) * xsq;
    }

  result = x * (xnum + a[3]) / (xden + b[3]);

  return result;
}

/*
 * Normal cdf for 0.66291 < fabs(x) < sqrt(32).
 */
static MpIeee gauss_medium(const MpIeee x)
{
  unsigned int  i;
  MpIeee temp=  MpIeee( "0.0" );
  MpIeee result=  MpIeee( "0.0" );
  MpIeee xnum;
  MpIeee xden;
  MpIeee absx;

  const MpIeee c[9] =  {
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8
  };
  const MpIeee d[8] =  {
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727
  };

  absx = fabs (x);

  xnum = c[8] * absx;
  xden = absx;

  for (i = 0; i < 7; i++)
    {
      xnum = (xnum + c[i]) * absx;
      xden = (xden + d[i]) * absx;
    }

  temp = (xnum + c[7]) / (xden + d[7]);

  result = get_del (x, temp);

  return result;
}

/*
 * Normal cdf for 
 * {sqrt(32) < x < GAUSS_XUPPER} union { GAUSS_XLOWER < x < -sqrt(32) }.
 */
static MpIeee gauss_large(const MpIeee x)
{
  int  i;
  MpIeee result;
  MpIeee xsq;
  MpIeee temp;
  MpIeee xnum;
  MpIeee xden;
  MpIeee absx;

  const MpIeee p[6] =  {
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303
  };
  const MpIeee q[5] =  {
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5
  };

  absx = fabs (x);
  xsq = MpIeee( "1.0" ) / (x * x);
  xnum = p[5] * xsq;
  xden = xsq;

  for (i = 0; i < 4; i++)
    {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }

  temp = xsq * (xnum + p[4]) / (xden + q[4]);
  temp = (M_1_SQRT2PI - temp) / absx;

  result = get_del (x, temp);

  return result;
}

MpIeee gsl_cdf_ugaussian_P(const MpIeee x)
{
  MpIeee result;
  MpIeee absx=  fabs (x);

  if (absx < GAUSS_EPSILON)
    {
      result = MpIeee( "0.5" );
      return result;
    }
  else if (absx < MpIeee( "0.66291" ))
    {
      result = MpIeee( "0.5" ) + gauss_small (x);
      return result;
    }
  else if (absx < SQRT32)
    {
      result = gauss_medium (x);

      if (x > 0.0)
        {
          result = MpIeee( "1.0" ) - result;
        }

      return result;
    }
  else if (x > GAUSS_XUPPER)
    {
      result = MpIeee( "1.0" );
      return result;
    }
  else if (x < GAUSS_XLOWER)
    {
      result = MpIeee( "0.0" );
      return result;
    }
  else
    {
      result = gauss_large (x);

      if (x > 0.0)
        {
          result = MpIeee( "1.0" ) - result;
        }
    }

  return result;
}

MpIeee gsl_cdf_ugaussian_Q(const MpIeee x)
{
  MpIeee result;
  MpIeee absx=  fabs (x);

  if (absx < GAUSS_EPSILON)
    {
      result = MpIeee( "0.5" );
      return result;
    }
  else if (absx < MpIeee( "0.66291" ))
    {
      result = gauss_small (x);

      if (x < 0.0)
        {
          result = fabs (result) + MpIeee( "0.5" );
        }
      else
        {
          result = MpIeee( "0.5" ) - result;
        }

      return result;
    }
  else if (absx < SQRT32)
    {
      result = gauss_medium (x);

      if (x < 0.0)
        {
          result = MpIeee( "1.0" ) - result;
        }

      return result;
    }
  else if (x > -(GAUSS_XLOWER))
    {
      result = MpIeee( "0.0" );
      return result;
    }
  else if (x < -(GAUSS_XUPPER))
    {
      result = MpIeee( "1.0" );
      return result;
    }
  else
    {
      result = gauss_large (x);

      if (x < 0.0)
        {
          result = MpIeee( "1.0" ) - result;
        }

    }

  return result;
}

MpIeee gsl_cdf_gaussian_P(const MpIeee x, const MpIeee sigma)
{
  return gsl_cdf_ugaussian_P (x / sigma);
}

MpIeee gsl_cdf_gaussian_Q(const MpIeee x, const MpIeee sigma)
{
  return gsl_cdf_ugaussian_Q (x / sigma);
}
