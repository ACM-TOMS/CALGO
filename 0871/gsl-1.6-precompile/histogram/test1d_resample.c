#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* histogram/test1d_resample.c
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
#include <math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

#include "urand.c"

void
test1d_resample (void)
{
  size_t i;
  int  status=  0;

  gsl_histogram *h;

  gsl_ieee_env_setup ();

  h = gsl_histogram_calloc_uniform (10, 0.0, 1.0);

  gsl_histogram_increment (h, 0.1);
  gsl_histogram_increment (h, 0.2);
  gsl_histogram_increment (h, 0.2);
  gsl_histogram_increment (h, 0.3);

  {
    gsl_histogram_pdf *p = gsl_histogram_pdf_alloc (10);

    gsl_histogram *hh = gsl_histogram_calloc_uniform (100, 0.0, 1.0);

    gsl_histogram_pdf_init (p, h);

    for (i = 0; i < 100000; i++)
      {
        MpIeee u=  urand();
        MpIeee x=  gsl_histogram_pdf_sample (p, u);
        gsl_histogram_increment (hh, x);
      }

    for (i = 0; i < 100; i++)
      {
        MpIeee y=  gsl_histogram_get (hh, i) / MpIeee( "2500" );
        MpIeee x;MpIeee  xmax;
        size_t k;
        MpIeee ya;

        gsl_histogram_get_range (hh, i, &x, &xmax);

        gsl_histogram_find (h, x, &k);
        ya = gsl_histogram_get (h, k);

        if (ya == MpIeee( "0" ))
          {
            if (y != MpIeee( "0" ))
              {
                {cout<<""<< (int)<<": "<<setiosflags((ios::floatfield))<< i;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" vs "<<setiosflags((ios::floatfield))<< y;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
                status = 1;
              }
          }
        else
          {
            MpIeee err=  MpIeee( "1" ) / sqrt (gsl_histogram_get (hh, i));
            MpIeee sigma=  fabs ((y - ya) / (ya * err));
            if (sigma > MpIeee( "3" ))
              {
                status = 1;
                {cout<<""<<setiosflags((ios::floatfield))<< y;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" vs "<<setiosflags((ios::floatfield))<< ya;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" err="<<setiosflags((ios::floatfield))<< err;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" sigma="<<setiosflags((ios::floatfield))<< sigma;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
              }
          }
      }

    gsl_histogram_pdf_free (p) ;
    gsl_histogram_free (hh);

    gsl_test (status, "gsl_histogram_pdf_sample within statistical errors");
  }

  gsl_histogram_free (h);
}
