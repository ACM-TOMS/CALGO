#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "demo_fn.h"
#include "demo_fn.c"

int
main (void)
{
  int  status;
  int  iter=  0;int   max_iter=  100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  MpIeee r=  MpIeee( "0" );MpIeee  r_expected=  sqrt (MpIeee( "5.0" ));
  MpIeee x_lo=  MpIeee( "0.0" );MpIeee  x_hi=  MpIeee( "5.0" );
  gsl_function F;
  struct quadratic_params params = {1.0, 0.0, -5.0};

  F.function = &quadratic;
  F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  {cout<<"using "<< 
          gsl_root_fsolver_name (s)<<" method\n";}

  {cout<<""<<setw(5)<<
          "iter";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" ["<<setw(9)<< "lower";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<", "<<setw(9)<< "upper";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<"] "<<setw(9)<< "root";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" "<<setw(10)<< 
          "err";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" "<<setw(9)<< "err(est)";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<"\n";}

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
                                       0, 0.001);

      if (status == GSL_SUCCESS)
        {cout<<"Converged:\n";}

      {cout<<""<<setw(5)<<
              iter;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" ["<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(7)<< x_lo;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<", "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(7)<< x_hi;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"] "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(7)<<
              r;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" "<<setiosflags((ios::showpos)|(ios::fixed & ios::floatfield))<<setprecision(7)<< r - r_expected;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::fixed & ios::floatfield))<<" "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(7)<< 
              x_hi - x_lo;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  return status;
}
