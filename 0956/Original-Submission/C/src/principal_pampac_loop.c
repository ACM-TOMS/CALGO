#include "mpi.h"
#include "pampac.h"

void
principal_pampac_loop (PTnode ** root_addr, options_struct * opts, int N_p) {
  double time_init, time_final, lambda, lambda_min, lambda_max, h, h_min;
  int lambda_index, global_iter, max_global_iter;
  char message[50];
  bool has_failed, has_completed;
  PTnode *root = *root_addr;
  debug_print (3, opts, __func__, "Beginning main computation...\n");

  /* Setting convenient aliases for optional parameters */
  lambda_min = opts->lambda_min;
  lambda_max = opts->lambda_max;
  lambda_index = opts->lambda_index;
  lambda = root->z[lambda_index];
  h_min = opts->h_min;
  h = root->h;
  max_global_iter = opts->max_global_iter;
  debug_print (0, opts, __func__, "h=%7.1e, h_min=%7.1e\n", h, h_min);
  debug_print (0, opts, __func__, "Maximum global iterations=%d\n",
               max_global_iter);

  global_iter = 0;
  has_failed = false;
  has_completed = (lambda <= lambda_min) || (lambda >= lambda_max);
  debug_print (0, opts, __func__, "Beginning global iteration.\n");

  /* Log image of initial node if necessary */
  visualize_tree (root, opts, "Initialization");

  time_init = MPI_Wtime ();
  while (!has_completed && !has_failed) {
    global_iter++;
    debug_print (5, opts, __func__, "while iteration %i:\n", global_iter);

    /* Spawn new nodes on tree at leaves if possible. */
    construct_predictor_nodes (root, opts);
    assign_processes (root, opts, N_p);
    assign_predictor_steps (root, opts);
    sprintf (message, "Iter %i: Assigned predictor steps.",
             global_iter);
    visualize_tree (root, opts, message);

    compute_corrector_steps (root, N_p);
    assess_residuals (root, opts);
    sprintf (message, "Iter %i: Corrector steps.", global_iter);
    visualize_tree (root, opts, message);

    prune_diverged_nodes (root, opts);
    sprintf (message, "Iter %i: Pruned DIVERGED nodes.", global_iter);
    visualize_tree (root, opts, message);

    construct_viable_paths (root, opts);
    choose_viable_paths (root, opts);
    sprintf (message, "Iter %i: Selected VIABLE paths.", global_iter);
    visualize_tree (root, opts, message);

    advance_root_node (&root, opts);
    sprintf (message, "Iter %i: Updated root node.", global_iter);
    visualize_tree (root, opts, message);

    lambda = root->z[lambda_index];
    has_completed = (lambda <= lambda_min) || (lambda >= lambda_max);
    h = root->h;

    has_failed = (global_iter >= max_global_iter) || (h < h_min);
    if (has_failed) {
      debug_print (0, opts, __func__, "Premature termination\n"); 
      debug_print (0, opts, __func__, 
                   "(global_iter = %i >= %i = max_global_iter) OR\n",
                   global_iter, max_global_iter);
      debug_print (0, opts, __func__,
                   "(h = %7.1e < %7.1e = h_min)\n", h, h_min);
      break;
    }
  }
  time_final = MPI_Wtime ();
  *root_addr = root;
  if (has_completed) {
    debug_print (0, opts, __func__, "Completed main loop\n");
    debug_print (0, opts, __func__, "%7.1e <= lambda = %7.1e <= %7.1e\n",
                 lambda_min, lambda, lambda_max); 
  }
  debug_print (0, opts, __func__, "Time elapsed = %g\n", time_final - time_init);
  debug_print (0, opts, __func__, "Global iterations = %d.\n", global_iter);
  return;
}
