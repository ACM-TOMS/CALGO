#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_rng.h>

gsl_rng * r;  /* global generator */

int
main (void)
{
  const gsl_rng_type * T;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  {cout<<"generator type: "<< gsl_rng_name (r)<<"\n";}
  {cout<<"seed = "<< gsl_rng_default_seed<<"\n";}
  {cout<<"first value = "<< gsl_rng_get (r)<<"\n";}
  return 0;
}
