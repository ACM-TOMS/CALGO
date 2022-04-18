#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_rng.h>

int
main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  int  i;int   n=  10;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (i = 0; i < n; i++) 
    {
      MpIeee u=  gsl_rng_uniform (r);
      {cout<<""<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(5)<< u;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
    }

  gsl_rng_free (r);

  return 0;
}
