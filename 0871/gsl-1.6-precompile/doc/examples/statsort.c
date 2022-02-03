#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

int
main(void)
{
  MpIeee data[5] =  {MpIeee( "17.2" ), MpIeee( "18.1" ), MpIeee( "16.5" ), MpIeee( "18.3" ), MpIeee( "12.6" )};
  MpIeee median;MpIeee  upperq;MpIeee  lowerq;

  {cout<<"Original dataset:  "<<setiosflags((ios::floatfield))<<
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

  gsl_sort (data, 1, 5);

  {cout<<"Sorted dataset: "<<setiosflags((ios::floatfield))<<
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

  median 
    = gsl_stats_median_from_sorted_data (data, 
                                         MpIeee( "1" ), MpIeee( "5" ));

  upperq 
    = gsl_stats_quantile_from_sorted_data (data, 
                                           MpIeee( "1" ), MpIeee( "5" ),
                                           MpIeee( "0.75" ));
  lowerq 
    = gsl_stats_quantile_from_sorted_data (data, 
                                           MpIeee( "1" ), MpIeee( "5" ),
                                           MpIeee( "0.25" ));

  {cout<<"The median is "<<setiosflags((ios::floatfield))<< median;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  {cout<<"The upper quartile is "<<setiosflags((ios::floatfield))<< upperq;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  {cout<<"The lower quartile is "<<setiosflags((ios::floatfield))<< lowerq;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  return 0;
}
