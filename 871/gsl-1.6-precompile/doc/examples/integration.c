#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

MpIeee f(MpIeee x, void * params) {
  MpIeee alpha=  *(MpIeee *) params;
  MpIeee f=  log(alpha*x) / sqrt(x);
  return f;
}

int
main (void)
{
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  
  MpIeee result;MpIeee  error;
  MpIeee expected=  -MpIeee( "4.0" );
  MpIeee alpha=  MpIeee( "1.0" );

  gsl_function F;
  F.function = &f;
  F.params = &alpha;

  gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
                        w, &result, &error); 

  {cout<<"result          = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(18)<< result;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"exact result    = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(18)<< expected;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"estimated error = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(18)<< error;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"actual error    = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(18)<< result - expected;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
  {cout<<"intervals =  "<< w->size<<"\n";}

  return 0;
}
