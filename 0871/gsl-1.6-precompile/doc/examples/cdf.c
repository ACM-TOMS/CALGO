#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_cdf.h>

int
main (void)
{
  MpIeee P;MpIeee  Q;
  MpIeee x=  MpIeee( "2.0" );

  P = gsl_cdf_ugaussian_P (x);
  {cout<<"prob(x < "<<setiosflags((ios::fixed & ios::floatfield))<< x;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<") = "<<setiosflags((ios::fixed & ios::floatfield))<< P;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}

  Q = gsl_cdf_ugaussian_Q (x);
  {cout<<"prob(x > "<<setiosflags((ios::fixed & ios::floatfield))<< x;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<") = "<<setiosflags((ios::fixed & ios::floatfield))<< Q;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}

  x = gsl_cdf_ugaussian_Pinv (P);
  {cout<<"Pinv("<<setiosflags((ios::fixed & ios::floatfield))<< P;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<") = "<<setiosflags((ios::fixed & ios::floatfield))<< x;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}

  x = gsl_cdf_ugaussian_Qinv (Q);
  {cout<<"Qinv("<<setiosflags((ios::fixed & ios::floatfield))<< Q;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<") = "<<setiosflags((ios::fixed & ios::floatfield))<< x;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}

  return 0;
}
