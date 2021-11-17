#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>

int
main (int argc, char **argv)
{
  int  i;int   n=  256;int   nc=  20;
  MpIeee *data=  malloc (n * sizeof (MpIeee));
  MpIeee *abscoeff=  malloc (n * sizeof (MpIeee));
  size_t *p = malloc (n * sizeof (size_t));

  FILE *f = fopen (argv[1], "r");
  for (i = 0; i < n; i++)
    {
      fscanf (f, "%lg", &data[i]);
    }
  fclose (f);

  {
    gsl_wavelet *w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
    gsl_wavelet_workspace *work = gsl_wavelet_workspace_alloc (n);

    gsl_wavelet_transform_forward (w, data, 1, n, work);

    for (i = 0; i < n; i++)
      {
        abscoeff[i] = fabs (data[i]);
      }

    gsl_sort_index (p, abscoeff, 1, n);

    for (i = 0; (i + nc) < n; i++)
      data[p[i]] = MpIeee( "0" );

    gsl_wavelet_transform_inverse (w, data, 1, n, work);
  }

  for (i = 0; i < n; i++)
    {
      {cout<<""<<setiosflags((ios::floatfield))<< data[i];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
    }
}
