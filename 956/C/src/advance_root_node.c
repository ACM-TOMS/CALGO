#include "pampac.h"
/**********************************************************************/
/* This function advances the root node in the event that it has only */
/* one child node and that node has converged (i.e., is CONVERGED).   */
/* This updating iterates to successive depths in the event that      */
/* there are several successive depths with solitary CONVERGED nodes. */
/**********************************************************************/
bool
advance_root_node (PTnode ** root_addr, options_struct * opts) {
  int k;
  bool has_solitary_child, has_converged, has_succeeded;
  PTnode *root_tmp = *root_addr;

  debug_print (3, opts, __func__, "Advancing root node of tree...\n");

  /*
   * The pointer root_tmp always points to the root node. The
   * associated root node is deleted iff it has a solitary CONVERGED 
   * child.
   */
  has_solitary_child = count_children (root_tmp) == 1;
  has_converged = root_tmp->state == CONVERGED;
  while ( has_solitary_child && has_converged ) {
    /* For loop to identify the index "k" of single non-NULL child. */
    for (k = 0; k < root_tmp->max_children; k++)
      if (root_tmp->child[k]!=NULL)
        break;
    /* Assertion to ensure child node is indeed CONVERGED */
    if (root_tmp->child[k]->state != CONVERGED)
      break;
    /* Dump z data from root node to disk before deleting. */
    has_succeeded = write_root_coordinates (root_tmp, opts);
    /* Update root_tmp's pointers to become new root node. */
    root_tmp = root_tmp->child[k];
    /* Carefully release pointers associated with old root node. */
    (*root_addr)->child[k] = NULL;
    free ((*root_addr)->z_init);
    (*root_addr)->z_init = NULL;
    free ((*root_addr)->T_init);
    (*root_addr)->T_init = NULL;
    free ((*root_addr)->child);
    (*root_addr)->child = NULL;
    /*
     * Reminder: we do not free *root_addr->z because the child's field
     * z_init points to this memory. Also, we release memory allocated
     * explicitly for the PTnode *root_addr itself.
     */
    free (*root_addr);
    *root_addr = root_tmp;
    debug_print (0, opts, __func__, "Advancing root node: lambda=%10g\n",
                 (*root_addr)->z[opts->lambda_index]);
    /* Check whether new root node has a solitary CONVERGED child. */
    has_solitary_child = count_children (root_tmp) == 1;
    has_converged = root_tmp->state == CONVERGED;
  }
  assign_depth (*root_addr, 0); /* Update to reflect changes in root node. */
  return has_succeeded;
}
