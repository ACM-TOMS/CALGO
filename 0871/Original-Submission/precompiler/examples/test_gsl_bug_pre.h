#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"


MpIeee test(){
  return MpIeee( "2" )
}

MpIeee *
gsl_vector_long_double_const_ptr(const gsl_vector_long_double * v, const size_t i)
{

/*
#if GSL_RANGE_CHECK
  if (i >= v->size)
    {
      GSL_ERROR_NULL ("index out of range", GSL_EINVAL);
    }
#endif
*/
  MpIeee b=  MpIeee( "2.0" );
  return  long MpIeee * (v->data + i * v->stride);
}


