#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>

int
main (void)
{
  MpIeee x=  MpIeee( "5.0" );
  MpIeee expected=  -MpIeee( "0.17759677131433830434739701" );
  
  MpIeee y=  gsl_sf_bessel_J0 (x);

  {cout<<"J0(5.0) = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(18)<< y;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"exact   = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(18)<< expected;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  return 0;
}
