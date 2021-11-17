#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_matrix.h>

int
main (void)
{
  int  i;int   j; 
  gsl_matrix * m = gsl_matrix_alloc (10, 3);
  
  for (i = 0; i < 10; i++)
    for (j = 0; j < 3; j++)
      gsl_matrix_set (m, i, j, 0.23 + 100*i + j);
  
  for (i = 0; i < 10; i++)
    for (j = 0; j < 3; j++)
      {cout<<"m("<< i<<","<< j<<") = "<<setiosflags((ios::floatfield))<< 
              gsl_matrix_get (m, i, j);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}

  return 0;
}
