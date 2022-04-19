#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* histogram/test2d_resample.c
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
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

#include "urand.c"

void
test2d_resample (void)
{
  size_t i, j;
  int  status=  0;
  MpIeee total=  MpIeee( "0" );
  size_t N = 200000;

  gsl_histogram2d *h;

  gsl_ieee_env_setup ();

  h = gsl_histogram2d_calloc_uniform (10, 10, 0.0, 1.0, 0.0, 1.0);

  for (i = 0; i < 10; i++)
    {
      for (j = 0; j < 10; j++)
        {
          MpIeee w=  MpIeee( "10.0" ) * i + j;
          total += w;
          gsl_histogram2d_accumulate (h, 0.1 * i, 0.1 * i, w);
        }
    }

  {
    gsl_histogram2d_pdf *p = gsl_histogram2d_pdf_alloc (10,10);

    gsl_histogram2d *hh = gsl_histogram2d_calloc_uniform (20, 20,
                                                          0.0, 1.0,
                                                          0.0, 1.0);

    gsl_histogram2d_pdf_init (p, h);

    for (i = 0; i < N; i++)
      {
        MpIeee u=  urand();
        MpIeee v=  urand();
        MpIeee x;MpIeee  y;
        status = gsl_histogram2d_pdf_sample (p, u, v, &x, &y);
        status = gsl_histogram2d_increment (hh, x, y);
      }

    status = 0;
    for (i = 0; i < 20; i++)
      {
        for (j = 0; j < 20; j++)
          {
            MpIeee z=  MpIeee( "4" ) * total * gsl_histogram2d_get (hh, i, j) / (MpIeee) N;
            size_t k1, k2;
            MpIeee ya;
            MpIeee x;MpIeee  xmax;MpIeee  y;MpIeee  ymax;

            gsl_histogram2d_get_xrange (hh, i, &x, &xmax);
            gsl_histogram2d_get_yrange (hh, j, &y, &ymax);

            gsl_histogram2d_find (h, x, y, &k1, &k2);
            ya = gsl_histogram2d_get (h, k1, k2);

            if (ya == MpIeee( "0" ))
              {
                if (z != MpIeee( "0" ))
                  {
                    status = 1;
                    {cout<<"("<< (int)<<","<<i<<"): "<<setiosflags((ios::floatfield))<< (int);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" vs "<<setiosflags((ios::floatfield))<<j;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
                  }
              }
            else
              {
                MpIeee err=  MpIeee( "1" ) / sqrt (gsl_histogram2d_get (hh, i, j));
                MpIeee sigma=  fabs ((z - ya) / (ya * err));
                if (sigma > MpIeee( "3" ))
                  {
                    status = 1;
                    {cout<<""<<setiosflags((ios::floatfield))<< z;
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
      }

    gsl_histogram2d_pdf_free (p) ;
    gsl_histogram2d_free (hh) ;
    
    gsl_test (status, "gsl_histogram2d_pdf_sample within statistical errors");
  }

  gsl_histogram2d_free (h) ;
}
