#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_linalg.h>

int
main (void)
{
  MpIeee a_data[] =  { MpIeee( "0.18" ), MpIeee( "0.60" ), MpIeee( "0.57" ), MpIeee( "0.96" ),
                      MpIeee( "0.41" ), MpIeee( "0.24" ), MpIeee( "0.99" ), MpIeee( "0.58" ),
                      MpIeee( "0.14" ), MpIeee( "0.30" ), MpIeee( "0.97" ), MpIeee( "0.66" ),
                      MpIeee( "0.51" ), MpIeee( "0.13" ), MpIeee( "0.19" ), MpIeee( "0.85" ) };

  MpIeee b_data[] =  { MpIeee( "1.0" ), MpIeee( "2.0" ), MpIeee( "3.0" ), MpIeee( "4.0" ) };

  gsl_matrix_view m 
    = gsl_matrix_view_array (a_data, 4, 4);

  gsl_vector_view b
    = gsl_vector_view_array (b_data, 4);

  gsl_vector *x = gsl_vector_alloc (4);
  
  int  s;

  gsl_permutation * p = gsl_permutation_alloc (4);

  gsl_linalg_LU_decomp (&m.matrix, p, &s);

  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

  {cout<<"x = \n";}
  gsl_vector_fprintf (stdout, x, "%g");

  gsl_permutation_free (p);
  return 0;
}
