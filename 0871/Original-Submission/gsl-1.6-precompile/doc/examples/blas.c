#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_blas.h>

int
main (void)
{
  MpIeee a[] =  { MpIeee( "0.11" ), MpIeee( "0.12" ), MpIeee( "0.13" ),
                 MpIeee( "0.21" ), MpIeee( "0.22" ), MpIeee( "0.23" ) };

  MpIeee b[] =  { MpIeee( "1011" ), MpIeee( "1012" ),
                 MpIeee( "1021" ), MpIeee( "1022" ),
                 MpIeee( "1031" ), MpIeee( "1032" ) };

  MpIeee c[] =  { MpIeee( "0.00" ), MpIeee( "0.00" ),
                 MpIeee( "0.00" ), MpIeee( "0.00" ) };

  gsl_matrix_view A = gsl_matrix_view_array(a, 2, 3);
  gsl_matrix_view B = gsl_matrix_view_array(b, 3, 2);
  gsl_matrix_view C = gsl_matrix_view_array(c, 2, 2);

  /* Compute C = A B */

  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &A.matrix, &B.matrix,
                  0.0, &C.matrix);

  {cout<<"[ "<<setiosflags((ios::floatfield))<< c[0];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", "<<setiosflags((ios::floatfield))<< c[1];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
  {cout<<"  "<<setiosflags((ios::floatfield))<< c[2];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<", "<<setiosflags((ios::floatfield))<< c[3];
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<" ]\n";}

  return 0;  
}
