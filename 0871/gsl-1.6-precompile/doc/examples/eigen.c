#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

int
main (void)
{
  MpIeee data[] =  { MpIeee( "1.0" )  , MpIeee( "1" )/MpIeee( "2.0" ), MpIeee( "1" )/MpIeee( "3.0" ), MpIeee( "1" )/MpIeee( "4.0" ),
                    MpIeee( "1" )/MpIeee( "2.0" ), MpIeee( "1" )/MpIeee( "3.0" ), MpIeee( "1" )/MpIeee( "4.0" ), MpIeee( "1" )/MpIeee( "5.0" ),
                    MpIeee( "1" )/MpIeee( "3.0" ), MpIeee( "1" )/MpIeee( "4.0" ), MpIeee( "1" )/MpIeee( "5.0" ), MpIeee( "1" )/MpIeee( "6.0" ),
                    MpIeee( "1" )/MpIeee( "4.0" ), MpIeee( "1" )/MpIeee( "5.0" ), MpIeee( "1" )/MpIeee( "6.0" ), MpIeee( "1" )/MpIeee( "7.0" ) };

  gsl_matrix_view m 
    = gsl_matrix_view_array (data, 4, 4);

  gsl_vector *eval = gsl_vector_alloc (4);
  gsl_matrix *evec = gsl_matrix_alloc (4, 4);

  gsl_eigen_symmv_workspace * w = 
    gsl_eigen_symmv_alloc (4);
  
  gsl_eigen_symmv (&m.matrix, eval, evec, w);

  gsl_eigen_symmv_free (w);

  gsl_eigen_symmv_sort (eval, evec, 
                        GSL_EIGEN_SORT_ABS_ASC);
  
  {
    int  i;

    for (i = 0; i < 4; i++)
      {
        MpIeee eval_i=  gsl_vector_get (eval, i);
        gsl_vector_view evec_i 
           = gsl_matrix_column (evec, i);

        {cout<<"eigenvalue = "<<setiosflags((ios::floatfield))<< eval_i;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
        {cout<<"eigenvector = \n";}
        gsl_vector_fprintf (stdout, 
                            &evec_i.vector, "%g");
      }
  }

  return 0;
}
