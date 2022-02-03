
double test(){
  return 2
}

long double *
gsl_vector_long_double_const_ptr (const gsl_vector_long_double * v, const size_t i)
{

/*
#if GSL_RANGE_CHECK
  if (i >= v->size)
    {
      GSL_ERROR_NULL ("index out of range", GSL_EINVAL);
    }
#endif
*/
  long double b = 2.0;
  return  long double * (v->data + i * v->stride);
}


