#include "pampac.h"
/**********************************************************************/
/* This routine functions primarily as a wrapper to the user's output */
/* routine write_coordinates. The primary purpose is to pass the data */
/* from the PTnode data structure to the user transparently.          */
/**********************************************************************/
bool
write_root_coordinates (PTnode *node, options_struct *opts) {
  debug_print (0, opts, __func__,
               "Writing point on curve to file, lambda = %8.1e\n",
               node->z[opts->lambda_index]);
  debug_print (5, opts, __func__,
               "node->N_dim=%i, node->z=%p\n", node->N_dim, node->z);
  debug_print (5, opts, __func__,
               "node->z[0]=%g, node->z[1]=%g...\n",
               node->z[0], node->z[1]);
  bool has_succeeded = write_coordinates (node->N_dim, node->z);
  return has_succeeded;
}
