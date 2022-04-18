#include <gsl_cblas.h>
#include "pampac.h"
/**********************************************************************/
/* Computes normalised secant vector from two points on a curve.      */
/**********************************************************************/
double
compute_secant_direction (PTnode *alpha, options_struct* opts) {
  int k, N;
  double *z, *z_old, *delta_z, delta_z_norm;
  debug_print (5, opts, __func__, "Computing secant direction.\n");
  N = alpha->N_dim;
  z = alpha->z;
  z_old = alpha->z_init;
  delta_z = alpha->T_init;
  for (k = 0; k < N; k++)
    delta_z[k] = z[k] - z_old[k];
  delta_z_norm = cblas_dnrm2 (N, delta_z, 1); /* BLAS routine for 2-norm */
  for (k = 0; k < N; k++)
    delta_z[k] /= delta_z_norm;
  return (delta_z_norm);
}
