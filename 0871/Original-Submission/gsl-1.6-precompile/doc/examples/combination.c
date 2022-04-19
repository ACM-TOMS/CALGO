#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_combination.h>

int 
main (void) 
{
  gsl_combination * c;
  size_t i;

  {cout<<"All subsets of {0,1,2,3} by size:\n";}
  for (i = 0; i <= 4; i++)
    {
      c = gsl_combination_calloc (4, i);
      do
        {
          {cout<<"{";}
          gsl_combination_fprintf (stdout, c, " %u");
          {cout<<" }\n";}
        }
      while (gsl_combination_next (c) == GSL_SUCCESS);
      gsl_combination_free (c);
    }

  return 0;
}
