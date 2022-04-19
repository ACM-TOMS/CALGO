#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

MpIeee f(MpIeee x, void * params)
{
  return pow (x, MpIeee( "1.5" ));
}

int
main (void)
{
  gsl_function F;
  MpIeee result;MpIeee  abserr;

  F.function = &f;
  F.params = 0;

  {cout<<"f(x) = x^(3/2)\n";}

  gsl_deriv_central (&F, 2.0, 1e-8, &result, &abserr);
  {cout<<"x = 2.0\n";}
  {cout<<"f'(x) = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(10)<< result;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" +/- "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(10)<< abserr;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"exact = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(10)<< 1.5 * sqrt(2.0);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n\n";}

  gsl_deriv_forward (&F, 0.0, 1e-8, &result, &abserr);
  {cout<<"x = 0.0\n";}
  {cout<<"f'(x) = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(10)<< result;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" +/- "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(10)<< abserr;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"exact = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(10)<< 0.0;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}

  return 0;
}
