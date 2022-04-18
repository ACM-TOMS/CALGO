#include <gsl_cblas.h>
#include "pampac.h"
bool initialize_secant (PTnode* root, options_struct *opts) {
  int count, k, N_dim, status;
  bool has_converged, has_failed;
  double *residual, r_nrm;

  debug_print (1, opts, __func__,
               "Computing secant direction using 2nd point.\n");

  /* Assume that node root has been initialized with a vector
   * root->z somehow. This function determines a nearby point on
   * the homotopy curve by perturbing the continuation parameter
   * a small amount and using corrector steps to return to the
   * curve. The two nearby points determine a secant direction. */
  N_dim = opts->N_dim;
  residual = malloc ((N_dim-1) * sizeof (double));
  /* Copy z into z_init */
  for (k=0; k<N_dim; k++)
    root->z_init[k] = root->z[k];

  /* Modify z_init by small perturbation in direction lambda_index */
  root->z_init[opts->lambda_index] -= (opts->lambda_dir) *
                                      (opts->delta_lambda);
  /* Initial secant direction is coordinate vector in the direction
   * of lambda_index (index of continuation parameter).  */
  for (k = 0; k < N_dim - 1; k++)
    root->T_init[k] = 0.e0;
  root->T_init[opts->lambda_index] = 1.0;

  debug_print (1, opts, __func__,
               "Iteration to get initial secant direction.\n");

  count = 0;
  has_converged = false;
  has_failed = false;
  while (!has_converged) {
    count++;
    has_failed = (count > opts->max_iter);
    if (has_failed) {
      debug_print (0, opts, __func__,
                   "Failed to determine initial secant direction.\n");
      debug_print (0, opts, __func__,
                   "Maximum of %d corrector iterations attained.\n",
                   opts->max_iter);
      debug_print (0, opts, __func__,
                   "Residual norm: %7.1e\n", opts->tol_residual);
      debug_print (0, opts, __func__,
                   "Desired residual norm: %7.1e\n",
                   opts->tol_residual);
      debug_print (0, opts, __func__, "Aborting processes.\n");
      return false;
    }
    status = single_corrector_step (N_dim, root->z_init, root->T_init);
    if (status!=0) {
      fprintf (stderr, "%s: Failed corrector step, status = %i\n",
               __func__, status);
      break;
    }
    status = compute_residual (N_dim, root->z_init, residual);
    if (status!=0) {
      fprintf (stderr, "%s: Error computing residual, status = %i\n",
               __func__, status);
      break;
    }
    r_nrm = cblas_dnrm2 (N_dim-1, residual, 1);
    debug_print (1, opts, __func__,
                 " count=%3d, residual norm=%7.1e.\n", count, r_nrm);
    has_converged = (r_nrm < opts->tol_residual);
  }
  free (residual);
  if (has_converged)
    compute_secant_direction (root, opts);
  else {
    debug_print (0, opts, __func__,
                 "Failed to determine initial secant direction\n");
    return false;
  }
  return true;
}
