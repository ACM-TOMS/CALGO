#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* histogram/pdf2d.c
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
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include "find.c"

int
 gsl_histogram2d_pdf_sample(const gsl_histogram2d_pdf * p,
                            MpIeee r1, MpIeee r2,
                            MpIeee *x, MpIeee *y)
{
  size_t k;
  int  status;

/* Wrap the exclusive top of the bin down to the inclusive bottom of
   the bin. Since this is a single point it should not affect the
   distribution. */

  if (r2 == MpIeee( "1.0" ))
    {
      r2 = MpIeee( "0.0" );
    }
  if (r1 == MpIeee( "1.0" ))
    {
      r1 = MpIeee( "0.0" );
    }

  status = find (p->nx * p->ny, p->sum, r1, &k);

  if (status)
    {
      GSL_ERROR ("cannot find r1 in cumulative pdf", GSL_EDOM);
    }
  else
    {
      size_t i = k / p->ny;
      size_t j = k - (i * p->ny);
      MpIeee delta=  (r1 - p->sum[k]) / (p->sum[k + 1] - p->sum[k]);
      *x = p->xrange[i] + delta * (p->xrange[i + 1] - p->xrange[i]);
      *y = p->yrange[j] + r2 * (p->yrange[j + 1] - p->yrange[j]);
      return GSL_SUCCESS;
    }
}

gsl_histogram2d_pdf *
gsl_histogram2d_pdf_alloc (const size_t nx, const size_t ny)
{
  const size_t n = nx * ny;

  gsl_histogram2d_pdf *p;

  if (n == 0)
    {
      GSL_ERROR_VAL ("histogram2d pdf length n must be positive integer",
                        GSL_EDOM, 0);
    }

  p = (gsl_histogram2d_pdf *) malloc (sizeof (gsl_histogram2d_pdf));

  if (p == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for histogram2d pdf struct",
                        GSL_ENOMEM, 0);
    }

  p->xrange = (MpIeee*) malloc ((nx + 1) * sizeof (MpIeee));

  if (p->xrange == 0)
    {
      free (p);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for histogram2d pdf xranges",
                        GSL_ENOMEM, 0);
    }

  p->yrange = (MpIeee*) malloc ((ny + 1) * sizeof (MpIeee));

  if (p->yrange == 0)
    {
      free (p->xrange);
      free (p);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for histogram2d pdf yranges",
                        GSL_ENOMEM, 0);
    }

  p->sum = (MpIeee*) malloc ((n + 1) * sizeof (MpIeee));

  if (p->sum == 0)
    {
      free (p->yrange);
      free (p->xrange);
      free (p);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for histogram2d pdf sums",
                        GSL_ENOMEM, 0);
    }

  p->nx = nx;
  p->ny = ny;

  return p;
}

int
 gsl_histogram2d_pdf_init(gsl_histogram2d_pdf * p, const gsl_histogram2d * h)
{
  size_t i;
  const size_t nx = p->nx;
  const size_t ny = p->ny;
  const size_t n = nx * ny;

  if (nx != h->nx || ny != h->ny)
    {
      GSL_ERROR ("histogram2d size must match pdf size", GSL_EDOM);
    }

  for (i = 0; i < n; i++)
    {
      if (h->bin[i] < 0)
        {
          GSL_ERROR ("histogram bins must be non-negative to compute"
                     "a probability distribution", GSL_EDOM);
        }
    }

  for (i = 0; i < nx + 1; i++)
    {
      p->xrange[i] = h->xrange[i];
    }

  for (i = 0; i < ny + 1; i++)
    {
      p->yrange[i] = h->yrange[i];
    }

  {
    MpIeee mean=  MpIeee( "0" );MpIeee  sum=  MpIeee( "0" );

    for (i = 0; i < n; i++)
      {
        mean += (h->bin[i] - mean) / ((MpIeee) (i + MpIeee( "1" )));
      }

    p->sum[0] = MpIeee( "0" );

    for (i = 0; i < n; i++)
      {
        sum += (h->bin[i] / mean) / n;
        p->sum[i + 1] = sum;
      }
  }

  return GSL_SUCCESS;
}


void
gsl_histogram2d_pdf_free (gsl_histogram2d_pdf * p)
{
  free (p->xrange);
  free (p->yrange);
  free (p->sum);
  free (p);
}
