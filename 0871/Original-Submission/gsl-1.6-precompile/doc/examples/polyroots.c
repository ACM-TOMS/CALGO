#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_poly.h>

int
main (void)
{
  int  i;
  /* coefficient of P(x) =  -1 + x^5  */
  MpIeee a[6] =  { -MpIeee( "1" ), MpIeee( "0" ), MpIeee( "0" ), MpIeee( "0" ), MpIeee( "0" ), MpIeee( "1" ) };  
  MpIeee z[10];

  gsl_poly_complex_workspace * w 
      = gsl_poly_complex_workspace_alloc (6);
  
  gsl_poly_complex_solve (a, 6, w, z);

  gsl_poly_complex_workspace_free (w);

  for (i = 0; i < 5; i++)
    {
      {cout<<"z"<< 
              i<<" = "<<setiosflags((ios::showpos)|(ios::fixed & ios::floatfield))<<setprecision(18)<< z[2*i];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::fixed & ios::floatfield))<<" "<<setiosflags((ios::showpos)|(ios::fixed & ios::floatfield))<<setprecision(18)<< z[2*i+1];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::showpos)|(ios::fixed & ios::floatfield))<<"\n";}
    }

  return 0;
}
