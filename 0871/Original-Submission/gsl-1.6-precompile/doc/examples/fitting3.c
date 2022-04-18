#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>

int
main (void)
{
  MpIeee x;
  const gsl_rng_type * T;
  gsl_rng * r;
  
  gsl_rng_env_setup ();
  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (x = MpIeee( "0.1" ); x < MpIeee( "2" ); x+= MpIeee( "0.1" ))
    {
      MpIeee y0=  exp (x);
      MpIeee sigma=  MpIeee( "0.1" ) * y0;
      MpIeee dy=  gsl_ran_gaussian (r, sigma);

      {cout<<""<<setiosflags((ios::floatfield))<< x;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<< y0 + dy;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<< sigma;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
    }
  return 0;
}
