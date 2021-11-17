#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

int
 func(MpIeee t, const MpIeee y[], MpIeee f[],
      void *params)
{
  MpIeee mu=  *(MpIeee *)params;
  f[0] = y[1];
  f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - MpIeee( "1" ));
  return GSL_SUCCESS;
}

int
 jac(MpIeee t, const MpIeee y[], MpIeee *dfdy, 
     MpIeee dfdt[], void *params)
{
  MpIeee mu=  *(MpIeee *)params;
  gsl_matrix_view dfdy_mat 
    = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix; 
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
  gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
  dfdt[0] = MpIeee( "0.0" );
  dfdt[1] = MpIeee( "0.0" );
  return GSL_SUCCESS;
}

int
main (void)
{
  const gsl_odeiv_step_type * T 
    = gsl_odeiv_step_rk8pd;

  gsl_odeiv_step * s 
    = gsl_odeiv_step_alloc (T, 2);
  gsl_odeiv_control * c 
    = gsl_odeiv_control_y_new (1e-6, 0.0);
  gsl_odeiv_evolve * e 
    = gsl_odeiv_evolve_alloc (2);

  MpIeee mu=  MpIeee( "10" );
  gsl_odeiv_system sys = {func, jac, 2, &mu};

  MpIeee t=  MpIeee( "0.0" );MpIeee  t1=  MpIeee( "100.0" );
  MpIeee h=  MpIeee( "1" )e-MpIeee( "6" );
  MpIeee y[2] =  { MpIeee( "1.0" ), MpIeee( "0.0" ) };

  while (t < t1)
    {
      int  status=  gsl_odeiv_evolve_apply (e, c, s,
                                           &sys, 
                                           &t, t1,
                                           &h, y);

      if (status != GSL_SUCCESS)
          break;

      {cout<<""<<setiosflags((ios::scientific & ios::floatfield))<<setprecision(5)<< t;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<" "<<setiosflags((ios::scientific & ios::floatfield))<<setprecision(5)<< y[0];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<" "<<setiosflags((ios::scientific & ios::floatfield))<<setprecision(5)<< y[1];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<"\n";}
    }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  return 0;
}
