#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_const_mksa.h>

int
main (void)
{
  MpIeee c=  GSL_CONST_MKSA_SPEED_OF_LIGHT;
  MpIeee au=  GSL_CONST_MKSA_ASTRONOMICAL_UNIT;
  MpIeee minutes=  GSL_CONST_MKSA_MINUTE;

  /* distance stored in meters */
  MpIeee r_earth=  MpIeee( "1.00" ) * au;  
  MpIeee r_mars=  MpIeee( "1.52" ) * au;

  MpIeee t_min;MpIeee  t_max;

  t_min = (r_mars - r_earth) / c;
  t_max = (r_mars + r_earth) / c;

  {cout<<"light travel time from Earth to Mars:\n";}
  {cout<<"minimum = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(1)<< t_min / minutes;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" minutes\n";}
  {cout<<"maximum = "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(1)<< t_max / minutes;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<" minutes\n";}

  return 0;
}
