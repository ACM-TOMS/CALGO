#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_statistics.h>

int
main(void)
{
  MpIeee data[5] =  {MpIeee( "17.2" ), MpIeee( "18.1" ), MpIeee( "16.5" ), MpIeee( "18.3" ), MpIeee( "12.6" )};
  MpIeee mean;MpIeee  variance;MpIeee  largest;MpIeee  smallest;

  mean     = gsl_stats_mean(data, MpIeee( "1" ), MpIeee( "5" ));
  variance = gsl_stats_variance(data, MpIeee( "1" ), MpIeee( "5" ));
  largest  = gsl_stats_max(data, MpIeee( "1" ), MpIeee( "5" ));
  smallest = gsl_stats_min(data, MpIeee( "1" ), MpIeee( "5" ));

  {cout<<"The dataset is "<<setiosflags((ios::floatfield))<<
         data[0];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", "<<setiosflags((ios::floatfield))<< data[1];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", "<<setiosflags((ios::floatfield))<< data[2];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", "<<setiosflags((ios::floatfield))<< data[3];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", "<<setiosflags((ios::floatfield))<< data[4];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}

  {cout<<"The sample mean is "<<setiosflags((ios::floatfield))<< mean;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  {cout<<"The estimated variance is "<<setiosflags((ios::floatfield))<< variance;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  {cout<<"The largest value is "<<setiosflags((ios::floatfield))<< largest;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  {cout<<"The smallest value is "<<setiosflags((ios::floatfield))<< smallest;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  return 0;
}
