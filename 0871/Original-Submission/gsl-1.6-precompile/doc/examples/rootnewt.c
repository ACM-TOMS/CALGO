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
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  MpIeee x0;MpIeee  x=  MpIeee( "5.0" );MpIeee  r_expected=  sqrt (MpIeee( "5.0" ));
  gsl_function_fdf FDF;
  struct quadratic_params params = {1.0, 0.0, -5.0};

  FDF.f = &quadratic;
  FDF.df = &quadratic_deriv;
  FDF.fdf = &quadratic_fdf;
  FDF.params = &params;

  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &FDF, x);

  {cout<<"using "<< 
          gsl_root_fdfsolver_name (s)<<" method\n";}

  {cout<<""<<setiosflags((ios::left))<<setw(5)<<
          "iter";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::left))<<" "<<setw(10)<< "root";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" "<<setw(10)<< "err";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" "<<setw(10)<< "err(est)";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<"\n";}
  do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, 1e-3);

      if (status == GSL_SUCCESS)
        {cout<<"Converged:\n";}

      {cout<<""<<setw(5)<<
              iter;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" "<<setiosflags((ios::fixed & ios::floatfield))<<setw(10)<<setprecision(7)<< x;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" "<<setiosflags((ios::showpos)|(ios::fixed & ios::floatfield))<<setw(10)<<setprecision(7)<< x - r_expected;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::fixed & ios::floatfield))<<" "<<setiosflags((ios::fixed & ios::floatfield))<<setw(10)<<setprecision(7)<< x - x0;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  return status;
}
