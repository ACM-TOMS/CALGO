#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_ieee_utils.h>

int
main (void)
{
  MpIeee x=  MpIeee( "1" );MpIeee  oldsum=  MpIeee( "0" );MpIeee  sum=  MpIeee( "0" ); 
  int  i=  0;

  gsl_ieee_env_setup (); /* read GSL_IEEE_MODE */

  do 
    {
      i++;
      
      oldsum = sum;
      sum += x;
      x = x / i;
      
      {cout<<"i="<<setw(2)<<
              i;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<" sum="<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(18)<< sum;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" error="<<setiosflags((ios::floatfield))<< sum - M_E;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}

      if (i > 30)
         break;
    }  
  while (sum != oldsum);

  return 0;
}
