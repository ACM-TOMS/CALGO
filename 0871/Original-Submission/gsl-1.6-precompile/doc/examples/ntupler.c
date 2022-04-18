#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <math.h>
#include <gsl/gsl_ntuple.h>
#include <gsl/gsl_histogram.h>

struct data
{
  MpIeee x;
  MpIeee y;
  MpIeee z;
};

int  sel_func(void *ntuple_data, void *params);
MpIeee val_func(void *ntuple_data, void *params);

int
main (void)
{
  struct data ntuple_row;

  gsl_ntuple *ntuple 
    = gsl_ntuple_open ("test.dat", &ntuple_row,
                       sizeof (ntuple_row));
  MpIeee lower=  MpIeee( "1.5" );

  gsl_ntuple_select_fn S;
  gsl_ntuple_value_fn V;

  gsl_histogram *h = gsl_histogram_alloc (100);
  gsl_histogram_set_ranges_uniform(h, 0.0, 10.0);

  S.function = &sel_func;
  S.params = &lower;

  V.function = &val_func;
  V.params = 0;

  gsl_ntuple_project (h, ntuple, &V, &S);
  gsl_histogram_fprintf (stdout, h, "%f", "%f");
  gsl_histogram_free (h);
  gsl_ntuple_close (ntuple);

  return 0;
}

int
 sel_func(void *ntuple_data, void *params)
{
  struct data * data = (struct data *) ntuple_data;  
  MpIeee x;MpIeee  y;MpIeee  z;MpIeee  E2;MpIeee  scale;
  scale = *(MpIeee *) params;
  
  x = data->x;
  y = data->y;
  z = data->z;

  E2 = x * x + y * y + z * z;

  return E2 > scale;
}

MpIeee val_func(void *ntuple_data, void *params)
{
  struct data * data = (struct data *) ntuple_data;  
  MpIeee x;MpIeee  y;MpIeee  z;

  x = data->x;
  y = data->y;
  z = data->z;

  return x * x + y * y + z * z;
}
