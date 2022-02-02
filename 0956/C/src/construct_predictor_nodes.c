#include "pampac.h"
/**********************************************************************/
/* This function traverses to the leaves of the tree and computes     */
/* predictor steps used to seed corrector iterations on other         */
/* processors. Nodes generated may exceed # of processes available.   */
/**********************************************************************/
void
construct_predictor_nodes (PTnode * alpha, options_struct * opts) {
  int k;
  bool is_max_depth, is_leaf_node, is_failed, has_data;
  PTnode *beta;

  /* Debug message prints only from root node */
  if (alpha->depth==0)
    debug_print (5, opts, __func__,
                 "Traversing tree, creating predictor nodes...\n");

  /* Predictor steps are constructed only when:
     - the node alpha is not FAILED
     - the node alpha's depth is below prescribed maximum
     - the node alpha is a leaf node
   */
  is_failed = (alpha->state == FAILED);
  /* Note: There are max_depth levels indexed from 0 */
  is_max_depth = (alpha->depth >= opts->max_depth-1);
  if (is_max_depth || is_failed)
    return;

  is_leaf_node = (count_children (alpha) == 0);
  if (is_leaf_node) {
    bool has_exceeded_max_step = false;
    /* Allocate memory for nodes for predictor steps & initialise. */
    for (k = 0; k < alpha->max_children; k++) {
      has_data = (alpha->T_init != NULL) && (alpha->z != NULL) &&
                 (alpha->z_init != NULL);
      double new_step = opts->scale_factors[k] * (alpha->h);
      bool is_too_long = (new_step > opts->h_max);
      /* Create at most one prediction with max step size. */
      if ((alpha->child[k] == NULL) && has_data &&
          !((has_exceeded_max_step) && (is_too_long))) {
        if (is_too_long) {
          new_step = opts->h_max;
          has_exceeded_max_step = true;
        }
        /* Allocate node and fill in scalar fields; assign actual
           process ID before allocating array fields of beta. */
        beta = initialize_PTnode (alpha->max_children);
        beta->N_dim = alpha->N_dim;
        beta->h_init = new_step;
        beta->h = new_step;
        beta->nu = 0;
        beta->nu_init = alpha->nu;
        beta->depth = alpha->depth + 1;
        alpha->child[k] = beta;
      }
    }
  } else {
    /* Loop over children, recursively descending to leaves. */
    for (k = 0; k < alpha->max_children; k++)
      if (alpha->child[k] != NULL)
        construct_predictor_nodes (alpha->child[k], opts);
  }
  return;
}
