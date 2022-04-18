#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>

int
main (void) 
{
  MpIeee f=  MpIeee( "1.0" )/MpIeee( "3.0" );
  MpIeee d=  MpIeee( "1.0" )/MpIeee( "3.0" );

  MpIeee fd=  f; /* promote from float to double */
  
  {cout<<" f=";} gsl_ieee_printf_float(&f); 
  {cout<<"\n";}

  {cout<<"fd=";} gsl_ieee_printf_double(&fd); 
  {cout<<"\n";}

  {cout<<" d=";} gsl_ieee_printf_double(&d); 
  {cout<<"\n";}

  return 0;
}
