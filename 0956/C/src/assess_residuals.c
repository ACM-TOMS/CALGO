#include <mpi.h>
#include "pampac.h"
/**********************************************************************/
/* This function traverses the tree, and, for each node in the tree,  */
/* receives the residual as computed by one of the processors. The    */
/* value received is used to assign states to the nodes which is the  */
/* basis of the decision-making algorithm for advancing on the curve. */
/**********************************************************************/
void
assess_residuals (PTnode * node, options_struct * opts) {
  int k, n_received;
  double residual_old, TOL, GAMMA, MU;
  bool has_failed, has_converged, has_almost_converged;
  MPI_Status status;

  /* Debug message prints only from root node */
  if (node->depth==0)
    debug_print (3, opts, __func__, 
                 "Traversing tree to assess state of each node...\n");

  TOL = opts->tol_residual;
  GAMMA = opts->gamma;
  MU = opts->mu;

  for (k = 0; k < node->max_children; k++)
    if (node->child[k] != NULL)
      assess_residuals (node->child[k], opts); /* Recursive call */

  if ((node->state==CONVERGED) || (node->pid<0))
    return;

  residual_old = node->res_norm; /* For measuring progress later */

  /* Corresponding MPI_Send calls in slave_process.
     Note: if node->pid==MPI_PROC_NULL, MPI_Recv will not hang. */

  MPI_Recv (node->z, node->N_dim, MPI_DOUBLE, node->pid,
            CONTINUE_TAG, MPI_COMM_WORLD, &status);
  /* Verifies that data was actually received */
  MPI_Get_count (&status, MPI_DOUBLE, &n_received);
  if (n_received==node->N_dim)
    node->nu++;
  else
    printf ("assess_residuals: n_received = %d != %d = node->N_dim\n",
            n_received, node->N_dim);
  MPI_Recv (&(node->res_norm), 1, MPI_DOUBLE, node->pid,
            CONTINUE_TAG, MPI_COMM_WORLD, &status);

  /* PROGRESSING    -> keep computing corrector steps (default)
     CONVERGED  -> converged; no more corrector steps needed
     CONVERGING -> almost converged; one more corrector step
     FAILED  -> insufficient progress or diverged         */

  has_failed = (node->nu > opts->max_iter) ||
               (log10(node->res_norm) > log10(MU) + log10(residual_old));
  has_converged = (node->res_norm < TOL);
  has_almost_converged = (GAMMA * log10 (node->res_norm) < log10 (TOL));
  node->state = PROGRESSING;
  if (has_converged)
    node->state = CONVERGED;
  else if (has_failed)
    node->state = FAILED;
  else if (has_almost_converged) {
    node->state = CONVERGING;
    debug_print (0, opts, __func__,
               "Residual=%7.1e/%-7.1e, nu=%d, h=%7.1e, state=",
               residual_old, node->res_norm, node->nu, node->h);
    print_state (node, stdout);
    printf("\n");
  }
  return;
}
