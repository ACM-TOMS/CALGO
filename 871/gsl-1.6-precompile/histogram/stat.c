#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* gsl_histogram_stat.c
 * Copyright (C) 2000  Simone Piccardi
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
 * File gsl_histogram_stat.c: 
 * Routines for statisticalcomputations on histograms. 
 * Need GSL library and header.
 * Contains the routines:
 * gsl_histogram_mean    compute histogram mean
 * gsl_histogram_sigma   compute histogram sigma
 *
 * Author: S. Piccardi
 * Jan. 2000
 *
 ***************************************************************/
#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>

/* FIXME: We skip negative values in the histogram h->bin[i] < 0,
   since those correspond to negative weights (BJG) */

MpIeee gsl_histogram_mean(const gsl_histogram * h)
{
  const size_t n = h->n;
  size_t i;

  /* Compute the bin-weighted arithmetic mean M of a histogram using the
     recurrence relation

     M(n) = M(n-1) + (x[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
     W(n) = W(n-1) + w(n)

   */

  MpIeee wmean=  MpIeee( "0" );
  MpIeee W=  MpIeee( "0" );

  for (i = 0; i < n; i++)
    {
      MpIeee xi=  (h->range[i + 1] + h->range[i]) / MpIeee( "2" );
      MpIeee wi=  h->bin[i];

      if (wi > MpIeee( "0" ))
        {
          W += wi;
          wmean += (xi - wmean) * (wi / W);
        }
    }

  return wmean;
}

MpIeee gsl_histogram_sigma(const gsl_histogram * h)
{
  const size_t n = h->n;
  size_t i;

  MpIeee wvariance=  MpIeee( "0" ) ;
  MpIeee wmean=  MpIeee( "0" );
  MpIeee W=  MpIeee( "0" );

  /* FIXME: should use a single pass formula here, as given in
     N.J.Higham 'Accuracy and Stability of Numerical Methods', p.12 */

  /* Compute the mean */

  for (i = 0; i < n; i++)
    {
      MpIeee xi=  (h->range[i + 1] + h->range[i]) / MpIeee( "2" );
      MpIeee wi=  h->bin[i];

      if (wi > MpIeee( "0" ))
        {
          W += wi;
          wmean += (xi - wmean) * (wi / W);
        }
    }

  /* Compute the variance */

  W = MpIeee( "0.0" );

  for (i = 0; i < n; i++)
    {
      MpIeee xi=  ((h->range[i + 1]) + (h->range[i])) / MpIeee( "2" );
      MpIeee wi=  h->bin[i];

      if (wi > MpIeee( "0" )) {
        const MpIeee delta=  (xi - wmean);
        W += wi ;
        wvariance += (delta * delta - wvariance) * (wi / W);
      }
    }

  {
    MpIeee sigma=  sqrt (wvariance) ;
    return sigma;
  }
}


/*
  sum up all bins of histogram
 */

MpIeee gsl_histogram_sum(const gsl_histogram * h)
{
  MpIeee sum= MpIeee( "0" );
  size_t i=0;
  size_t n;
  n=h->n;

  while(i < n)
    sum += h->bin[i++];

  return sum;
}

