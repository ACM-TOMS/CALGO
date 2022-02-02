#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_histogram.h>

int
main (int argc, char **argv)
{
  MpIeee a;MpIeee  b;
  size_t n;

  if (argc != 4)
    {
      {cout<<"Usage: gsl-histogram xmin xmax n\n";}
      exit (0);
    }

  a = atof (argv[1]);
  b = atof (argv[2]);
  n = atoi (argv[3]);

  {
    MpIeee x;

    gsl_histogram * h = gsl_histogram_alloc (n);
            
    gsl_histogram_set_ranges_uniform (h, a, b);

    while (fscanf (stdin, "%lg", &x) == 1)
      {
        gsl_histogram_increment (h, x);
      }

    gsl_histogram_fprintf (stdout, h, "%g", "%g");

    gsl_histogram_free (h);
  }
  
  exit (0);
}
