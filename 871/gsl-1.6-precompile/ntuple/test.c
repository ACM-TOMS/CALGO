#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_ntuple.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

struct data
{
  int  num;
  MpIeee x;
  MpIeee y;
  MpIeee z;
};
int  sel_func(void *ntuple_data, void * params);
MpIeee val_func(void *ntuple_data, void * params);

int
main (void)
{
  struct data ntuple_row;
  int  i;

  MpIeee x[1000];MpIeee  y[1000];MpIeee  z[1000];MpIeee  f[100];

  gsl_ntuple_select_fn S;
  gsl_ntuple_value_fn V;
  
  MpIeee scale=  MpIeee( "1.5" );
  
  gsl_ieee_env_setup ();

  S.function = &sel_func;
  S.params = &scale;
  
  V.function = &val_func;
  V.params = &scale;

  {
    gsl_ntuple *ntuple = gsl_ntuple_create ("test.dat", &ntuple_row, 
                                            sizeof (ntuple_row));

    int  status=  0;

    for (i = 0; i < 100; i++) f[i] = MpIeee( "0" );
    
    for (i = 0; i < 1000; i++)
      {
        MpIeee xi=  MpIeee( "1.0" ) / (i + MpIeee( "1.5" ));
        MpIeee yi=  xi * xi ;
        MpIeee zi=  xi * xi * xi;
        
        ntuple_row.x = xi;
        ntuple_row.y = yi;
        ntuple_row.z = zi;
        ntuple_row.num = i;
        
        x[i] = xi; y[i] = yi; z[i] = zi;
        
        if (xi * scale < MpIeee( "0.1" ))
          {
            MpIeee v=  xi + yi + zi;
            int  k=  (int)(100.0*v*scale);
            f[k]++;
          }

        /* printf ("x,y,z = %f,%f,%f; n=%x \n", ntuple_row.x,
           ntuple_row.y, ntuple_row.z, ntuple_row.num); */
        
        {
          int  s=  gsl_ntuple_bookdata (ntuple);

          if (s != GSL_SUCCESS)
            {
              status = 1;
            }
        }
      }
    
    gsl_ntuple_close (ntuple);

    gsl_test (status, "writing ntuples");
  }

  {
    gsl_ntuple *ntuple = gsl_ntuple_open ("test.dat", &ntuple_row, 
                                          sizeof (ntuple_row));
    int  status=  0;

    for (i = 0; i < 1000; i++)
      {
        gsl_ntuple_read (ntuple);

        status = (ntuple_row.num != i);
        status |= (ntuple_row.x != x[i]);
        status |= (ntuple_row.y != y[i]);
        status |= (ntuple_row.z != z[i]);

        /* printf ("x,y,z = %f,%f,%f; n=%d\n", ntuple_row.x,
                ntuple_row.y, ntuple_row.z, ntuple_row.num); */
      }
    gsl_ntuple_close (ntuple);

    gsl_test (status, "reading ntuples");
  }    

  {
    int  status=  0;

    gsl_ntuple *ntuple = gsl_ntuple_open ("test.dat", &ntuple_row, 
                                          sizeof (ntuple_row));

    gsl_histogram *h = gsl_histogram_calloc_uniform (100, 0., 1.);

    gsl_ntuple_project (h, ntuple, &V, &S);

    gsl_ntuple_close (ntuple);

    /* gsl_histogram_fprintf (stdout, h, "%f", "%f"); */

    for (i = 0; i < 100; i++)
      {
        /* printf ("h  %g f  %g\n", h->bin[i], f[i]); */

        if (h->bin[i] != f[i])
          {
            status = 1;
          }
      }

    gsl_test (status, "histogramming ntuples");

    gsl_histogram_free (h);
  }

  exit (gsl_test_summary());
}

int
 sel_func(void *ntuple_data, void * params)
{
  MpIeee x;MpIeee  y;MpIeee  z;MpIeee  scale;
  scale = *(MpIeee *)params;

  x = ((struct data *) ntuple_data)->x;
  y = ((struct data *) ntuple_data)->y;
  z = ((struct data *) ntuple_data)->z;

  return (x*scale < 0.1);
}

MpIeee val_func(void *ntuple_data, void * params)
{
  MpIeee x;MpIeee  y;MpIeee  z;MpIeee  scale;
  scale = *(MpIeee *)params;

  x = ((struct data *) ntuple_data)->x;
  y = ((struct data *) ntuple_data)->y;
  z = ((struct data *) ntuple_data)->z;

  return (x + y + z) * scale;
}
