#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
#include <gsl/gsl_vector.h>

int
main (void)
{
  int  i; 
  gsl_vector * v = gsl_vector_alloc (10);

  {  
     FILE * f = fopen ("test.dat", "r");
     gsl_vector_fscanf (f, v);
     fclose (f);
  }

  for (i = 0; i < 10; i++)
    {
      {cout<<""<<setiosflags((ios::floatfield))<< gsl_vector_get(v, i);
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::floatfield))<<"\n";}
    }

  return 0;
}
