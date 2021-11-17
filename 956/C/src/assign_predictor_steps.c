#include "pampac.h"
/**********************************************************************/
/* This function traverses carries out a breadth-first traversal of   */
/* the tree to assign data for predictor steps at leaf nodes where    */
/* the nodes require initialisation.                                  */
/**********************************************************************/
void
assign_predictor_steps (PTnode *root, options_struct *opts) {
  Queue q;
  PTnode *alpha, *beta;
  bool is_leaf, has_done_work, has_pid;
  int j, k, N_dim;

  N_dim = root->N_dim;
  /* Breadth-first traversal of tree */
  init_queue (&q);
  enqueue (&q, root);
  while (!empty_queue (&q)) {
    alpha = front_of_queue (&q);
    dequeue (&q); /* Mark alpha as visited by dequeuing */
    /* Examine children of alpha: append non-leaf nodes to end of
    queue, process leaf nodes */
    for (k=0; k<alpha->max_children; k++) {
      beta = alpha->child[k];
      if (beta != NULL) {
        is_leaf = (count_children(beta) == 0);
        if (!is_leaf)
          enqueue (&q, beta);
        else {
          /* has_done_work prevents erasing useful data. */
          has_done_work = (beta->nu>0);
          has_pid = (beta->pid>0);
          if (!has_done_work && has_pid) {
            if (beta->z==NULL)
              beta->z = malloc (N_dim * sizeof (double));
            /* z_init: shallow copy of z from parent */
            if (beta->z_init==NULL)
              beta->z_init = alpha->z;
            if (beta->T_init==NULL)
              beta->T_init = malloc (N_dim * sizeof (double));
            compute_secant_direction (alpha, opts);
            for (j=0; j<N_dim; j++) {
              beta->T_init[j] = alpha->T_init[j];
              beta->z[j] = beta->z_init[j] +
                           beta->h * beta->T_init[j];
            }
          }
        }
      }
    }
  }
  return;
}
