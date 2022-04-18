#include "pampac.h"
/**********************************************************************/
/* This routine simply traverses a tree of PTnodes, printing contents */
/* of each node to stdout as it goes.                                 */
/**********************************************************************/
void
print_tree (PTnode *alpha) {
  int k;
  if (alpha->child != NULL)
    for (k = 0; k < alpha->max_children; k++)
      if (alpha->child[k] != NULL)
        print_tree (alpha->child[k]);
  printf ("\n");
  print_PTnode (alpha);
  return;
}
