#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

int
main (void)
{
  MpIeee x=  MpIeee( "5.0" );
  gsl_sf_result result;

  MpIeee expected=  -MpIeee( "0.17759677131433830434739701" );
  
  int  status=  gsl_sf_bessel_J0_e (x, &result);

  {cout<<"status  = "<< gsl_strerror(status)<<"\n";}
  {cout<<"J0(5.0) = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(18)<<
          "      +/- % .18f\n";
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"exact   = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(18)<< expected;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  return status;
}
