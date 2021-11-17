#include "pampac.h"
bool create_root_node (PTnode** root_addr, options_struct *opts) {
  PTnode *root;
  int N_dim = opts->N_dim;
  bool is_allocated;

  root = initialize_PTnode (opts->max_children); /* Allocate root node */
  is_allocated = (root!=NULL);
  if (!is_allocated) {
    *root_addr = NULL;
    return false;
  }
  /* Fill in fields of root node appropriately from options */
  root->N_dim = N_dim;
  root->depth = 0;
  root->state = CONVERGED;
  root->nu = 0;
  root->h_init = 0.; /* The intial root has no parent */
  root->h = opts->h_init;
  /* Attempt to allocate memory for fields z, z_init, & T_init. */
  root->z = malloc (N_dim * sizeof (double));
  is_allocated = is_allocated && (root->z!=NULL);
  root->z_init = malloc (N_dim * sizeof (double));
  is_allocated = is_allocated && (root->z_init!=NULL);
  root->T_init = malloc (N_dim * sizeof (double));
  is_allocated = is_allocated && (root->T_init!=NULL);
  if (!is_allocated) {
    *root_addr = NULL;
    return false;
  }
  /* Pointer root will be deallocated on return to calling stack frame.
   * As such, its value must be copied explicitly to *root_addr to
   * ensure that the space allocated is preserved. */
  *root_addr = root;
  if (opts->verbose>4)
    print_PTnode (root);
  return true;
}
