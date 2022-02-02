#include "mpi.h"
#include "pampac.h"
#include <gsl_cblas.h>
/**********************************************************************/
/* Main routine executed by master processor. The bulk of the work is */
/* in the routine principal_pampac_loop; the remaining code is        */
/* largely for preparation and clean-up.                              */
/**********************************************************************/
void
master_process (int N_p, options_struct *opts) {
  PTnode *root = NULL;
  bool has_succeeded;

  /* Verification: Allocate root_node successfully */
  has_succeeded = create_root_node (&root, opts);
  if (!has_succeeded) {
    debug_print (0, opts, __func__,
                 "Failed to allocate memory for root node.\n");
    debug_print (0, opts, __func__,
                 "Terminating...\n");
    goto cleanup;
  }

  /* Verification: Load initial point from disk */
  has_succeeded = load_initial_coordinates (root, opts);
  if (!has_succeeded) {
    debug_print (0, opts, __func__,
                 "Failed to read first point into root node.\n");
    debug_print (0, opts, __func__,
                 "Terminating...\n");
    goto cleanup;
  }

  if (opts->verbose>0) {
    double res_nrm;
    int lambda_index = opts->lambda_index;
    debug_print (1, opts, __func__,
                 "Loaded first point from file.\n");
    int status = compute_residual (root->N_dim, root->z, root->T_init);
    if (status!=0) {
		fprintf (stderr, "%s: Error in computing residual.\n", __func__);
		fprintf (stderr, "%s: error status = %i\n", __func__, status);
		goto cleanup;
	}
    res_nrm = cblas_dnrm2 (root->N_dim-1, root->T_init, 1);
    debug_print (1, opts, __func__,
                 "Initial residual = %7.1e\n", res_nrm);
    debug_print (1, opts, __func__,
                 "Initial lambda = %8.1e\n", root->z[lambda_index]);
  }

  /* Verification: Determine initial secant direction */
  has_succeeded = initialize_secant (root, opts);
  if (!has_succeeded) {
    debug_print (0, opts, __func__,
                 "Failed to determine initial secant direction.\n");
    debug_print (0, opts, __func__,
                 "Maximum of %d iterations exceeded.\n", opts->max_iter);
    debug_print (0, opts, __func__,
                 "Desired residual tolerance: %7.1e\n",opts->tol_residual);
    debug_print (0, opts, __func__,
                 "Terminating.\n");
    goto cleanup;
  }
  principal_pampac_loop (&root, opts, N_p);
  write_root_coordinates (root, opts);

cleanup:
  debug_print (0, opts, __func__, "Shutting down slave processes.\n");
  stop_slaves (N_p);
  debug_print (0, opts, __func__, "Cleaning up memory.\n");

  delete_options (opts);

  if (root!=NULL)
    delete_tree (root, opts);
  return;
}
