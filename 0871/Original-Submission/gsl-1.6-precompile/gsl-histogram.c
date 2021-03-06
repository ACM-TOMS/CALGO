#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* histogram/gsl-histogram.c
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
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_histogram.h>

int
main (int argc, char **argv)
{
  MpIeee a=  MpIeee( "0.0" );MpIeee  b=  MpIeee( "1.0" );
  size_t n = 10;

  if (argc != 3 && argc !=4)
    {
      {cout<<"Usage: gsl-histogram xmin xmax [n]\n";}
      {cout<<"Computes a histogram of the data on stdin using n bins from xmin to xmax.\n";}
      exit (0);
    }

  a = atof (argv[1]);
  b = atof (argv[2]);

  if (argc == 4) 
    {
      n = atoi (argv[3]);
    }
  else
    {
      n = (int)(b - a) ;
    }

  {
    MpIeee x;
    gsl_histogram *h = gsl_histogram_alloc (n);

    gsl_histogram_set_ranges_uniform (h, a, b);

    while (fscanf(stdin, "%lg", &x) == 1)
      {
        gsl_histogram_increment(h, x);
      }

#ifdef DISPLAY_STATS
    {
      MpIeee mean=  gsl_histogram_mean (h);
      MpIeee sigma=  gsl_histogram_sigma (h);
      fprintf (stdout, "# mean = %g\n", mean);
      fprintf (stdout, "# sigma = %g\n", sigma);
    }
#endif

    gsl_histogram_fprintf (stdout, h, "%g", "%g");

    gsl_histogram_free (h);
  }

  return 0;
}
