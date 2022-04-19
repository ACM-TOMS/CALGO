#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_fit.h>

int
main (void)
{
  int  i;int   n=  4;
  MpIeee x[4] =  { MpIeee( "1970" ), MpIeee( "1980" ), MpIeee( "1990" ), MpIeee( "2000" ) };
  MpIeee y[4] =  {   MpIeee( "12" ),   MpIeee( "11" ),   MpIeee( "14" ),   MpIeee( "13" ) };
  MpIeee w[4] =  {  MpIeee( "0.1" ),  MpIeee( "0.2" ),  MpIeee( "0.3" ),  MpIeee( "0.4" ) };

  MpIeee c0;MpIeee  c1;MpIeee  cov00;MpIeee  cov01;MpIeee  cov11;MpIeee  chisq;

  gsl_fit_wlinear (x, 1, w, 1, y, 1, n, 
                   &c0, &c1, &cov00, &cov01, &cov11, 
                   &chisq);

  {cout<<"# best fit: Y = "<<setiosflags((ios::floatfield))<< c0;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" + "<<setiosflags((ios::floatfield))<< c1;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" X\n";}
  {cout<<"# covariance matrix:\n";}
  {cout<<"# [ "<<setiosflags((ios::floatfield))<< 
          cov00;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", "<<setiosflags((ios::floatfield))<< cov01;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n#   "<<setiosflags((ios::floatfield))<< cov01;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", "<<setiosflags((ios::floatfield))<< cov11;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"]\n";}
  {cout<<"# chisq = "<<setiosflags((ios::floatfield))<< chisq;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}

  for (i = 0; i < n; i++)
    {cout<<"data: "<<setiosflags((ios::floatfield))<< 
                   x[i];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<< y[i];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<< 1/sqrt(w[i]);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}

  {cout<<"\n";}

  for (i = -30; i < 130; i++)
    {
      MpIeee xf=  x[0] + (i/MpIeee( "100.0" )) * (x[n-1] - x[0]);
      MpIeee yf;MpIeee  yf_err;

      gsl_fit_linear_est (xf, 
                          c0, c1, 
                          cov00, cov01, cov11, 
                          &yf, &yf_err);

      {cout<<"fit: "<<setiosflags((ios::floatfield))<< xf;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<< yf;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
      {cout<<"hi : "<<setiosflags((ios::floatfield))<< xf;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<< yf + yf_err;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
      {cout<<"lo : "<<setiosflags((ios::floatfield))<< xf;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<< yf - yf_err;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
    }
  return 0;
}
