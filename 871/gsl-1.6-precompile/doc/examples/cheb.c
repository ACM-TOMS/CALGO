#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>

MpIeee f(MpIeee x, void *p)
{
  if (x < MpIeee( "0.5" ))
    return MpIeee( "0.25" );
  else
    return MpIeee( "0.75" );
}

int
main (void)
{
  int  i;int   n=  10000; 

  gsl_cheb_series *cs = gsl_cheb_alloc (40);

  gsl_function F;

  F.function = f;
  F.params = 0;

  gsl_cheb_init (cs, &F, 0.0, 1.0);

  for (i = 0; i < n; i++)
    {
      MpIeee x=  i / (MpIeee)n;
      MpIeee r10=  gsl_cheb_eval_n (cs, MpIeee( "10" ), x);
      MpIeee r40=  gsl_cheb_eval (cs, x);
      {cout<<""<<setiosflags((ios::floatfield))<< 
              x;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<< GSL_FN_EVAL (&F, x);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<<, r10;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" "<<setiosflags((ios::floatfield))<< r40;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
    }

  gsl_cheb_free (cs);

  return 0;
}
