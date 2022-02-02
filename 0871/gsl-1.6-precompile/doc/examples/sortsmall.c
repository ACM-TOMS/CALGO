#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>

int
main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  size_t i, k = 5, N = 100000;

  MpIeee * x=  malloc (N * sizeof(MpIeee));
  MpIeee * small=  malloc (k * sizeof(MpIeee));

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (i = 0; i < N; i++)
    {
      x[i] = gsl_rng_uniform(r);
    }

  gsl_sort_smallest (small, k, x, 1, N);

  {cout<<""<< k<<" smallest values from "<< N<<"\n";}

  for (i = 0; i < k; i++)
    {
      {cout<<""<< i<<": "<<setiosflags((ios::fixed & ios::floatfield))<<setprecision(18)<< small[i];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"\n";}
    }
  return 0;
}
