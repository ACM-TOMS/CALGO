#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

MpIeee fn1(MpIeee x, void * params)
{
  return cos(x) + MpIeee( "1.0" );
}

int
main (void)
{
  int  status;
  int  iter=  0;int   max_iter=  100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  MpIeee m=  MpIeee( "2.0" );MpIeee  m_expected=  M_PI;
  MpIeee a=  MpIeee( "0.0" );MpIeee  b=  MpIeee( "6.0" );
  gsl_function F;

  F.function = &fn1;
  F.params = 0;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m, a, b);

  {cout<<"using "<<
          gsl_min_fminimizer_name (s)<<" method\n";}

  {cout<<""<<setw(5)<<
          "iter";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" ["<<setw(9)<< "lower";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<", "<<setw(9)<< "upper";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<"] "<<setw(9)<< "min";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" "<<setw(10)<<
          "err";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" "<<setw(9)<< "err(est)";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<"\n";}

  {cout<<""<<setw(5)<<
          iter;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" ["<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(7)<< a;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<", "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(7)<< b;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"] "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(7)<<
          m;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" "<<setiosflags((ios::showpos)|(ios::fixed & ios::floatfield))<<setprecision(7)<< m - m_expected;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::fixed & ios::floatfield))<<" "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(7)<< b - a;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}

  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status 
        = gsl_min_test_interval (a, b, 0.001, 0.0);

      if (status == GSL_SUCCESS)
        {cout<<"Converged:\n";}

      {cout<<""<<setw(5)<<
              "%.7f %.7f %+.7f %.7f\n";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" ["<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(7)<<
              iter;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<", "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(7)<< a;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"] ";}
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  return status;
}
