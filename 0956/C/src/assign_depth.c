#include "pampac.h"
/**********************************************************************/
/* Traverse tree, labelling depth of nodes en route.                  */
/* As nodes are deleted and children replace their parents, the depth */
/* of nodes on the tree changes; hence, it is necessary to call this  */
/* function (starting from the root node) to label the depth of each  */
/* node in the tree correctly after updating the root node.           */
/**********************************************************************/
void
assign_depth (PTnode * node, int depth) {
  int k;
  for (k = 0; k < node->max_children; k++) {
    if (node->child[k] != NULL) {
      assign_depth (node->child[k], depth + 1);
    }
  }
  node->depth = depth;
}
