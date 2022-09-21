#include "cals_blas.h"

void set_threads(int threads) {
  omp_set_num_threads(threads);
#if CALS_MKL
  mkl_set_num_threads(threads);
#elif CALS_BLIS
  bli_thread_set_num_threads(threads);
#elif CALS_OPENBLAS
  openblas_set_num_threads(threads);
#elif CALS_MATLAB
  omp_set_num_threads(threads);
#endif
}

int get_threads() {
#if CALS_MKL
  return mkl_get_max_threads();
#elif CALS_BLIS
  return bli_thread_get_num_threads();
#elif CALS_OPENBLAS
  return openblas_get_num_threads();
#elif CALS_MATLAB
  return omp_get_max_threads();
#endif
}
