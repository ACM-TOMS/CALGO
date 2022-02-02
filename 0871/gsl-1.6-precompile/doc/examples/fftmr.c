#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

int
main (void)
{
  int  i;
  const int n = 630;
  MpIeee data[2*n];

  gsl_fft_complex_wavetable * wavetable;
  gsl_fft_complex_workspace * workspace;

  for (i = 0; i < n; i++)
    {
      REAL(data,i) = 0.0;
      IMAG(data,i) = 0.0;
    }

  data[0] = MpIeee( "1.0" );

  for (i = 1; i <= 10; i++)
    {
      REAL(data,i) = REAL(data,n-i) = 1.0;
    }

  for (i = 0; i < n; i++)
    {
      {cout<<""<< i<<": "<<setiosflags((ios::scientific & ios::floatfield))<< REAL(data,i);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<" "<<setiosflags((ios::scientific & ios::floatfield))<<, 
                                IMAG(data,i);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<"\n";}
    }
  {cout<<"\n";}

  wavetable = gsl_fft_complex_wavetable_alloc (n);
  workspace = gsl_fft_complex_workspace_alloc (n);

  for (i = 0; i < wavetable->nf; i++)
    {
       {cout<<"# factor "<< i<<": "<< 
               wavetable->factor[i]<<"\n";}
    }

  gsl_fft_complex_forward (data, 1, n, 
                           wavetable, workspace);

  for (i = 0; i < n; i++)
    {
      {cout<<""<< i<<": "<<setiosflags((ios::scientific & ios::floatfield))<< REAL(data,i);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<" "<<setiosflags((ios::scientific & ios::floatfield))<<, 
                                IMAG(data,i);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<"\n";}
    }

  gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);
  return 0;
}
