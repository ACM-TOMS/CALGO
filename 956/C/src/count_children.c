#include "pampac.h"
/**********************************************************************/
/* Accumulates total number of non-NULL children of a tree node.      */
/**********************************************************************/
int
count_children (PTnode * alpha) {
  int k, count = 0;
  if (alpha->child != NULL)
    for (k = 0; k < alpha->max_children; k++) {
      if (alpha->child[k] != NULL)
        count++;
    }
  return count;
}
