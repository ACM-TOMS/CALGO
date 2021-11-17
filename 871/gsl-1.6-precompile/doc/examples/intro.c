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
  MpIeee y=  gsl_sf_bessel_J0 (x);
  {cout<<"J0("<<setiosflags((ios::floatfield))<< x;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<") = "<<setiosflags((ios::scientific & ios::floatfield))<<setprecision(18)<< y;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::scientific & ios::floatfield))<<"\n";}
  return 0;
}
