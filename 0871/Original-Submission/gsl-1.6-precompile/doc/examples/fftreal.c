#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

int
main (void)
{
  int  i;int   n=  100;
  MpIeee data[n];

  gsl_fft_real_wavetable * real;
  gsl_fft_halfcomplex_wavetable * hc;
  gsl_fft_real_workspace * work;

  for (i = 0; i < n; i++)
    {
      data[i] = MpIeee( "0.0" );
    }

  for (i = n / 3; i < 2 * n / 3; i++)
    {
      data[i] = MpIeee( "1.0" );
    }

  for (i = 0; i < n; i++)
    {
      {cout<<""<< i<<": "<<setiosflags((ios::scientific & ios::floatfield))<< data[i];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<"\n";}
    }
  {cout<<"\n";}

  work = gsl_fft_real_workspace_alloc (n);
  real = gsl_fft_real_wavetable_alloc (n);

  gsl_fft_real_transform (data, 1, n, 
                          real, work);

  gsl_fft_real_wavetable_free (real);

  for (i = 11; i < n; i++)
    {
      data[i] = MpIeee( "0" );
    }

  hc = gsl_fft_halfcomplex_wavetable_alloc (n);

  gsl_fft_halfcomplex_inverse (data, 1, n, 
                               hc, work);
  gsl_fft_halfcomplex_wavetable_free (hc);

  for (i = 0; i < n; i++)
    {
      {cout<<""<< i<<": "<<setiosflags((ios::scientific & ios::floatfield))<< data[i];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<"\n";}
    }

  gsl_fft_real_workspace_free (work);
  return 0;
}
