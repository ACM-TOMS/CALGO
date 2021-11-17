#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* histogram/stat2d.c
 * Copyright (C) 2002  Achim Gaedke
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/***************************************************************
 *
 * File histogram/stat2d.c:
 * Routine to return statistical values of the content of a 2D hisogram. 
 *
 * Contains the routines:
 * gsl_histogram2d_sum sum up all bin values
 * gsl_histogram2d_xmean determine mean of x values
 * gsl_histogram2d_ymean determine mean of y values
 *
 * Author: Achim Gaedke Achim.Gaedke@zpr.uni-koeln.de
 * Jan. 2002
 *
 ***************************************************************/

#include <config.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram2d.h>

/*
  sum up all bins of histogram2d
 */

MpIeee gsl_histogram2d_sum(const gsl_histogram2d * h)
{
  const size_t n = h->nx * h->ny;
  MpIeee sum=  MpIeee( "0" );
  size_t i = 0;

  while (i < n)
    sum += h->bin[i++];

  return sum;
}

MpIeee gsl_histogram2d_xmean(const gsl_histogram2d * h)
{
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  size_t i;
  size_t j;

  /* Compute the bin-weighted arithmetic mean M of a histogram using the
     recurrence relation

     M(n) = M(n-1) + (x[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
     W(n) = W(n-1) + w(n)

   */

  MpIeee wmean=  MpIeee( "0" );
  MpIeee W=  MpIeee( "0" );

  for (i = 0; i < nx; i++)
    {
      MpIeee xi=  (h->xrange[i + 1] + h->xrange[i]) / MpIeee( "2.0" );
      MpIeee wi=  MpIeee( "0" );

      for (j = 0; j < ny; j++)
        {
          MpIeee wij=  h->bin[i * ny + j];
          if (wij > MpIeee( "0" ))
            wi += wij;
        }
      if (wi > MpIeee( "0" ))
        {
          W += wi;
          wmean += (xi - wmean) * (wi / W);
        }
    }

  return wmean;
}

MpIeee gsl_histogram2d_ymean(const gsl_histogram2d * h)
{
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  size_t i;
  size_t j;

  /* Compute the bin-weighted arithmetic mean M of a histogram using the
     recurrence relation

     M(n) = M(n-1) + (x[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
     W(n) = W(n-1) + w(n)

   */

  MpIeee wmean=  MpIeee( "0" );
  MpIeee W=  MpIeee( "0" );

  for (j = 0; j < ny; j++)
    {
      MpIeee yj=  (h->yrange[j + 1] + h->yrange[j]) / MpIeee( "2.0" );
      MpIeee wj=  MpIeee( "0" );

      for (i = 0; i < nx; i++)
        {
          MpIeee wij=  h->bin[i * ny + j];
          if (wij > MpIeee( "0" ))
            wj += wij;
        }

      if (wj > MpIeee( "0" ))
        {
          W += wj;
          wmean += (yj - wmean) * (wj / W);
        }
    }

  return wmean;
}

MpIeee gsl_histogram2d_xsigma(const gsl_histogram2d * h)
{
  const MpIeee xmean=  gsl_histogram2d_xmean (h);
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  size_t i;
  size_t j;

  /* Compute the bin-weighted arithmetic mean M of a histogram using the
     recurrence relation

     M(n) = M(n-1) + (x[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
     W(n) = W(n-1) + w(n)

   */

  MpIeee wvariance=  MpIeee( "0" );
  MpIeee W=  MpIeee( "0" );

  for (i = 0; i < nx; i++)
    {
      MpIeee xi=  (h->xrange[i + 1] + h->xrange[i]) / MpIeee( "2" ) - xmean;
      MpIeee wi=  MpIeee( "0" );

      for (j = 0; j < ny; j++)
        {
          MpIeee wij=  h->bin[i * ny + j];
          if (wij > MpIeee( "0" ))
            wi += wij;
        }

      if (wi > MpIeee( "0" ))
        {
          W += wi;
          wvariance += ((xi * xi) - wvariance) * (wi / W);
        }
    }

  {
    MpIeee xsigma=  sqrt (wvariance);
    return xsigma;
  }
}

MpIeee gsl_histogram2d_ysigma(const gsl_histogram2d * h)
{
  const MpIeee ymean=  gsl_histogram2d_ymean (h);
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  size_t i;
  size_t j;

  /* Compute the bin-weighted arithmetic mean M of a histogram using the
     recurrence relation

     M(n) = M(n-1) + (x[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
     W(n) = W(n-1) + w(n)

   */

  MpIeee wvariance=  MpIeee( "0" );
  MpIeee W=  MpIeee( "0" );

  for (j = 0; j < ny; j++)
    {
      MpIeee yj=  (h->yrange[j + 1] + h->yrange[j]) / MpIeee( "2.0" ) - ymean;
      MpIeee wj=  MpIeee( "0" );

      for (i = 0; i < nx; i++)
        {
          MpIeee wij=  h->bin[i * ny + j];
          if (wij > MpIeee( "0" ))
            wj += wij;
        }
      if (wj > MpIeee( "0" ))
        {
          W += wj;
          wvariance += ((yj * yj) - wvariance) * (wj / W);
        }
    }

  {
    MpIeee ysigma=  sqrt (wvariance);
    return ysigma;
  }
}

MpIeee gsl_histogram2d_cov(const gsl_histogram2d * h)
{
  const MpIeee xmean=  gsl_histogram2d_xmean (h);
  const MpIeee ymean=  gsl_histogram2d_ymean (h);
  const size_t nx = h->nx;
  const size_t ny = h->ny;
  size_t i;
  size_t j;

  /* Compute the bin-weighted arithmetic mean M of a histogram using the
     recurrence relation

     M(n) = M(n-1) + (x[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
     W(n) = W(n-1) + w(n)

   */

  MpIeee wcovariance=  MpIeee( "0" );
  MpIeee W=  MpIeee( "0" );

  for (j = 0; j < ny; j++)
    {
      for (i = 0; i < nx; i++)
        {
          MpIeee xi=  (h->xrange[i + 1] + h->xrange[i]) / MpIeee( "2.0" ) - xmean;
          MpIeee yj=  (h->yrange[j + 1] + h->yrange[j]) / MpIeee( "2.0" ) - ymean;
          MpIeee wij=  h->bin[i * ny + j];

          if (wij > MpIeee( "0" ))
            {
              W += wij;
              wcovariance += ((xi * yj) - wcovariance) * (wij / W);
            }
        }
    }

  return wcovariance;

}
