#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* fft/compare_source.c
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

#include "compare.h"

int
 FUNCTION(compare_complex,results) (const char *name_a, const BASE a[],
                                   const char *name_b, const BASE b[],
                                   size_t stride, size_t n,
                                   const MpIeee allowed_ticks)
{
  size_t i;
  MpIeee ticks;MpIeee  max_ticks=  MpIeee( "0" );
  MpIeee dr;MpIeee  di;
  const char *flag;

  for (i = 0; i < n; i++)
    {
      dr = b[2*stride*i] - a[2*stride*i];
      di = b[2*stride*i+1] - a[2*stride*i+1];
      ticks = (fabs (dr) + fabs (di)) / BASE_EPSILON;
      if (ticks > max_ticks)
        {
          max_ticks = ticks;
        }
    }

  if (max_ticks < allowed_ticks)
    {
      return 0;
    }

  {cout<<"\n"<< name_a<<" vs "<< name_b<<" : max_ticks = "<<setiosflags((ios::fixed & ios::floatfield))<< max_ticks;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}

  for (i = 0; i < n; i++)
    {
      dr = b[2*stride*i] - a[2*stride*i];
      di = b[2*stride*i+1] - a[2*stride*i+1];
      ticks = (fabs (dr) + fabs (di)) / BASE_EPSILON;

      if (ticks > MpIeee( "1000" ))
        {
          flag = "***";
        }
      else
        {
          flag = "";
        }

      {cout<<""<<setw(15)<< name_a;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<": "<< (int)<<"  "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(16)<<i;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(16)<<
              a[2*stride*i];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" "<< a[2*stride*i+1]<<"\n";}
      {cout<<""<<setw(15)<< name_b;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<": "<< (int)<<"  "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(16)<<i;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(16)<<
              b[2*stride*i];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" "<<setiosflags((ios::scientific & ios::floatfield))<< b[2*stride*i+1];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<" "<< ticks<<"\n";}
    }

  return -1;
}


int
 FUNCTION(compare_real,results) (const char *name_a, const BASE a[],
                                const char *name_b, const BASE b[],
                                size_t stride, size_t n,
                                const MpIeee allowed_ticks)
{
  size_t i;
  MpIeee ticks;MpIeee  max_ticks=  MpIeee( "0" );
  MpIeee dr;
  const char *flag;

  for (i = 0; i < n; i++)
    {
      dr = b[stride*i] - a[stride*i];
      ticks = fabs (dr) / BASE_EPSILON;
      if (ticks > max_ticks)
        {
          max_ticks = ticks;
        }
    }

  if (max_ticks < allowed_ticks)
    {
      return 0;
    }

  {cout<<"\n"<< name_a<<" vs "<< name_b<<" : max_ticks = "<<setiosflags((ios::fixed & ios::floatfield))<< max_ticks;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}

  for (i = 0; i < n; i++)
    {
      dr = b[stride*i] - a[stride*i];
      ticks = fabs (dr) / BASE_EPSILON;

      if (ticks > MpIeee( "1000" ))
        {
          flag = "***";
        }
      else
        {
          flag = "";
        }

      {cout<<""<<setw(15)<< name_a;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<": "<< (int)<<"  "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(16)<<i;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" "<< 
              a[stride*i]<<"\n";}
      {cout<<""<<setw(15)<< name_b;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<": "<< (int)<<"  "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(16)<<i;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" "<<setiosflags((ios::scientific & ios::floatfield))<< 
              b[stride*i];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<" "<< ticks<<"\n";}
    }

  return -1;
}
