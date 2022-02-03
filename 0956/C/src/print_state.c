#include "pampac.h"
/**********************************************************************/
/* Given a PTnode, print its state to stdout.                        */
/**********************************************************************/
void
print_state (PTnode *node, FILE *file) {
  switch (node->state) {
  case FAILED:
    fprintf (file, "FAILED");
    break;
  case CONVERGED:
    fprintf (file, "CONVERGED");
    break;
  case CONVERGING:
    fprintf (file, "CONVERGING");
    break;
  case PROGRESSING:
    fprintf (file, "PROGRESSING");
    break;
  default:
    break;
  }
}
