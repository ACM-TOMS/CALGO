#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

int
main (void)
{
  const gsl_odeiv_step_type * T 
    = gsl_odeiv_step_rk4;

  gsl_odeiv_step * s 
    = gsl_odeiv_step_alloc (T, 2);

  MpIeee mu=  MpIeee( "10" );
  gsl_odeiv_system sys = {func, jac, 2, &mu};

  MpIeee t=  MpIeee( "0.0" );MpIeee  t1=  MpIeee( "100.0" );
  MpIeee h=  MpIeee( "1" )e-MpIeee( "2" );
  MpIeee y[2] =  { MpIeee( "1.0" ), MpIeee( "0.0" ) };MpIeee  y_err[2];
  MpIeee dydt_in[2];MpIeee  dydt_out[2];

  /* initialise dydt_in from system parameters */
  GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);

  while (t < t1)
    {
      int  status=  gsl_odeiv_step_apply (s, t, h, 
                                         y, y_err, 
                                         dydt_in, 
                                         dydt_out, 
                                         &sys);

      if (status != GSL_SUCCESS)
          break;

      dydt_in[0] = dydt_out[0];
      dydt_in[1] = dydt_out[1];

      t += h;

      {cout<<""<<setiosflags((ios::scientific & ios::floatfield))<<setprecision(5)<< t;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<" "<<setiosflags((ios::scientific & ios::floatfield))<<setprecision(5)<< y[0];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<" "<<setiosflags((ios::scientific & ios::floatfield))<<setprecision(5)<< y[1];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<"\n";}
    }

  gsl_odeiv_step_free (s);
  return 0;
}
