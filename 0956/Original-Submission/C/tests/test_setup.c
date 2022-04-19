#include "test_pampac.h"

/* A sample tree to use for testing pampac library routines. */

options_struct
set_options(int depth, int max_children) {
  options_struct opts;
  int k;
  opts.N_dim = 6;
  opts.lambda_min = 0.0;
  opts.lambda_max = 1.0;
  opts.lambda_index = opts.N_dim-1;
  opts.lambda_dir = 1;
  opts.delta_lambda = 1.0e-5;
  opts.h_min = 1.0e-10;
  opts.h_max = 1.0e2;
  opts.h_init = 1.0e-1;
  opts.max_iter = 5;
  opts.tol_residual = 1.0e-6;
  opts.mu = 0.5;
  opts.gamma = 2.0;
  opts.max_depth = depth;
  opts.max_children = max_children;
  opts.scale_factors = malloc (max_children * sizeof (double));
  for (k=0; k<opts.max_children; k++)
    opts.scale_factors[k] = 0.5 + k*0.25;
  opts.input_filename = "input.txt";
  opts.tree_base_filename = "./data/trees";
  opts.tree_filename_num = 0;
  opts.verbose = 0;
  return (opts);
}

void free_options (options_struct *opts) {
  if (opts->scale_factors!=NULL)
    free (opts->scale_factors);
  return;
}

PTnode* make_tree (options_struct *opts) {
  int k;
  PTnode *root = initialize_PTnode (opts->max_children);
  root->N_dim = opts->N_dim;
  root->state = CONVERGED;
  root->h = opts->h_init;
  root->nu = 0;
  root->max_children = opts->max_children;
  root->depth = 0;
  root->z = malloc (root->N_dim * sizeof (double));
  root->T_init = malloc (root->N_dim * sizeof (double));
  root->z_init = malloc (root->N_dim * sizeof (double));
  for (k=0; k<root->N_dim; k++) {
    root->z[k] = (double)(k+2);
    root->z_init[k] = (double)(k+1);
  }
  compute_secant_direction (root, opts);
  /* Generate new levels of tree (nodes only) */
  for (k=1; k<opts->max_depth; k++) {
    construct_predictor_nodes (root, opts);
    assign_processes (root, opts, 50);  /* Assume large number of processors */
    assign_predictor_steps (root, opts);
  }
  return (root);
}


void tree43(PTnode* root, options_struct *opts) {
  /* Given tree of depth 4 with 3 children per node, assigns various
     states to the nodes and prunes a few subtrees. */
  root->state = CONVERGED;
  root->nu = 4;

  root->child[0]->state = CONVERGED;
  root->child[1]->state = CONVERGED;
  root->child[2]->state = CONVERGED;

  root->child[0]->nu = 4;
  root->child[1]->nu = 5;
  root->child[2]->nu = 3;

  root->child[0]->child[0]->state = CONVERGED;
  root->child[0]->child[1]->state = FAILED;
  root->child[0]->child[2]->state = CONVERGING;

  root->child[0]->child[0]->nu = 4;
  root->child[0]->child[1]->nu = 3;
  root->child[0]->child[2]->nu = 3;

  root->child[1]->child[0]->state = FAILED;
  delete_tree (root->child[1]->child[1], opts);
  root->child[1]->child[1]=NULL;
  root->child[1]->child[2]->state = FAILED;

  root->child[1]->child[0]->nu = 4;
  root->child[1]->child[2]->nu = 3;

  root->child[2]->child[0]->state = CONVERGING;
  root->child[2]->child[1]->state = PROGRESSING;
  root->child[2]->child[2]->state = FAILED;

  root->child[2]->child[0]->nu = 4;
  root->child[2]->child[1]->nu = 2;
  root->child[2]->child[2]->nu = 3;

  root->child[0]->child[0]->child[0]->state = FAILED;
  delete_tree (root->child[0]->child[0]->child[1], opts);
  root->child[0]->child[0]->child[1]=NULL;
  root->child[0]->child[0]->child[2]->state = FAILED;

  root->child[0]->child[0]->child[0]->nu = 3;
  root->child[0]->child[0]->child[2]->nu = 3;

  root->child[0]->child[1]->child[0]->state = PROGRESSING;
  root->child[0]->child[1]->child[1]->state = CONVERGING;
  root->child[0]->child[1]->child[2]->state = PROGRESSING;

  root->child[0]->child[1]->child[0]->nu = 3;
  root->child[0]->child[1]->child[1]->nu = 3;
  root->child[0]->child[1]->child[2]->nu = 3;

  root->child[0]->child[2]->child[0]->state = PROGRESSING;
  root->child[0]->child[2]->child[1]->state = FAILED;
  root->child[0]->child[2]->child[2]->state = CONVERGING;

  root->child[0]->child[2]->child[0]->nu = 2;
  root->child[0]->child[2]->child[1]->nu = 2;
  root->child[0]->child[2]->child[2]->nu = 2;

  root->child[1]->child[0]->child[0]->state = PROGRESSING;
  root->child[1]->child[0]->child[1]->state = CONVERGING;
  root->child[1]->child[0]->child[2]->state = PROGRESSING;

  root->child[1]->child[0]->child[0]->nu = 2;
  root->child[1]->child[0]->child[1]->nu = 4;
  root->child[1]->child[0]->child[2]->nu = 2;

  root->child[1]->child[2]->child[0]->state = PROGRESSING;
  root->child[1]->child[2]->child[1]->state = CONVERGED;
  root->child[1]->child[2]->child[2]->state = CONVERGING;

  root->child[1]->child[2]->child[0]->nu = 3;
  root->child[1]->child[2]->child[1]->nu = 5;
  root->child[1]->child[2]->child[2]->nu = 2;

  root->child[2]->child[0]->child[0]->state = CONVERGED;
  root->child[2]->child[0]->child[1]->state = CONVERGING;
  root->child[2]->child[0]->child[2]->state = PROGRESSING;

  root->child[2]->child[0]->child[0]->nu = 5;
  root->child[2]->child[0]->child[1]->nu = 2;
  root->child[2]->child[0]->child[2]->nu = 2;

  root->child[2]->child[1]->child[0]->state = FAILED;
  root->child[2]->child[1]->child[1]->state = CONVERGED;
  delete_tree (root->child[2]->child[1]->child[2], opts);
  root->child[2]->child[1]->child[2] = NULL;

  root->child[0]->child[2]->child[0]->nu = 2;
  root->child[0]->child[2]->child[1]->nu = 5;

  root->child[2]->child[2]->child[0]->state = CONVERGING;
  root->child[2]->child[2]->child[1]->state = CONVERGED;
  root->child[2]->child[2]->child[2]->state = FAILED;

  root->child[2]->child[2]->child[0]->nu = 4;
  root->child[2]->child[2]->child[0]->nu = 4;
  root->child[2]->child[2]->child[0]->nu = 2;

  return;
}
