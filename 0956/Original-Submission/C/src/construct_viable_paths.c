#include "pampac.h"
/**********************************************************************/
/* This function recursively traverses the tree, computing lengths of */
/* valid (containing all CONVERGED nodes) and viable (containing either   */
/* CONVERGED or CONVERGING nodes) paths to the leaves of the tree.            */
/**********************************************************************/
void
construct_viable_paths (PTnode *alpha, options_struct* opts) {
  int k, valid_index = -1, viable_index = -1;
  bool is_leaf, is_valid, is_viable, is_child_valid, is_child_viable;
  double valid_path_length = 0.0, viable_path_length = 0.0, length;
  PTnode *beta;

  is_leaf = (count_children(alpha)==0);
  is_valid = (alpha->state==CONVERGED);
  is_viable = is_valid || (alpha->state==CONVERGING);

  /* Default values NaN/-1; algorithm overwrites as necessary. */
  alpha->valid_path_length = NAN;
  alpha->valid_index = -1;
  alpha->nu_valid = -1;
  alpha->viable_path_length = NAN;
  alpha->viable_index = -1;
  alpha->nu_viable = -1;

  if (is_valid) {
    alpha->valid_path_length = alpha->h_init;
    alpha->nu_valid = alpha->nu;
  } /* In case alpha is the end of a valid path and
       is not a leaf node. */
  if (is_viable) {
    alpha->viable_path_length = alpha->h_init;
    alpha->nu_viable = alpha->nu;
  } /* In case alpha is the end of a viable path and is not a leaf node. */

  /* Leaf nodes: Zero path lengths as appropriate. */
  if (!is_leaf) {
    /* For non-leaf nodes, first recurse to leaf nodes, updating
    path lengths en route. Second, for each valid/viable child,
    add step-length to valid/viable path lengths provided alpha
    is also valid/viable. */
    for (k=0; k<alpha->max_children; k++) {
      beta = alpha->child[k];
      if (beta!=NULL) {
        construct_viable_paths (beta, opts);  /* recursive call */
        is_child_valid = (beta->state==CONVERGED);
        is_child_viable = is_child_valid ||
                          (beta->state==CONVERGING);
        if (is_viable && is_child_viable) {
          length = beta->viable_path_length + alpha->h_init;
          /* Test below fails if either side is NaN. */
          if (viable_path_length < length) {
            viable_path_length = length;
            viable_index = k;
          }
        }
        if (is_valid && is_child_valid) {
          length = beta->valid_path_length + alpha->h_init;
          /* Test below fails if either side is NaN. */
          if (valid_path_length < length) {
            valid_path_length = length;
            valid_index = k;
          }
        }
      }
    }
  }

  /* By default, valid/viable path lengths are NaN (i.e., when no paths
     include node alpha). If the locally computed valid/viable path
     length including alpha is strictly positive, update accordingly. */
  if (viable_path_length>0.0) {
    alpha->viable_path_length = viable_path_length;
    alpha->viable_index = viable_index;
    beta = alpha->child[viable_index];
    alpha->nu_viable = max(alpha->nu, beta->nu_init + beta->nu_viable);
  }
  if (valid_path_length>0.0) {
    alpha->valid_path_length = valid_path_length;
    alpha->valid_index = valid_index;
    beta = alpha->child[valid_index];
    alpha->nu_valid = max(alpha->nu, beta->nu_init + beta->nu_valid);
  }

  return;
}
