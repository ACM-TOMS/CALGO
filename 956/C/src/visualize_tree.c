#include <string.h>
#include "pampac.h"
/**********************************************************************/
/* This routine traverses a tree of PTnodes, generating output in the */
/* dot language (compatible with the GraphViz family of graph tools). */
/**********************************************************************/
void
visualize_tree (PTnode *root, options_struct *opts, const char* title) {
  int k, n_chars;
  Queue q;
  PTnode *alpha, *beta;
  char *filename;
  FILE *out_file;

  if (opts->verbose <=2)
    return;

  /* Parameter tree_filename_num is negative when a previously
   * attempted file-write failed (in which case, don't do more). */
  if (opts->tree_filename_num<0)
    return;

  /* Extract location of base output filename from opts */
  n_chars = strlen (opts->tree_base_filename) + 5 + 4;
  filename = malloc (n_chars * sizeof (char));
  /* Note: n_chars accounts for 5 chars for numeric field & 4 chars
   * for the additional characters "-" and ".gv". If the "%05d" is
   * modified in the next line, n_chars should be modified also.   */
  sprintf (filename, "%s-%05d.gv",
           opts->tree_base_filename, opts->tree_filename_num);
  out_file = fopen (filename, "w");

  free (filename); /* Must release before test that follows... */

  if (out_file == NULL) {
    debug_print (0, opts, __func__,
                 "Warning, error opening %s.\n", filename);
    /* Set signal to abort future calls to visualize_tree */
    opts->tree_filename_num = -1;
    return;
  }

  /* Start writing information to dot file */
  fprintf (out_file, "digraph G {\n");
  fprintf (out_file,
           "node [shape=circle,style=filled,fontcolor=black];\n");
  /* Breadth-first traversal of tree */
  init_queue (&q);
  enqueue (&q, root);
  while (!empty_queue (&q)) {
    alpha = front_of_queue (&q);
    dequeue (&q); /* Mark alpha as visited by dequeuing. */

    /* Print node information about alpha into GraphViz file. */
    fprintf (out_file, "%d [label=\"%d\\nnu=%d\\nPID=%d\", color=",
             alpha->label, alpha->label, alpha->nu, alpha->pid);
    switch (alpha->state) {
    case FAILED:
      fprintf (out_file, "BLACK");
      break;
    case CONVERGED:
      fprintf (out_file, "GREEN");
      break;
    case CONVERGING:
      fprintf (out_file, "YELLOW");
      break;
    case PROGRESSING:
      fprintf (out_file, "RED");
      break;
    default:
      break;
    }
    fprintf (out_file, "];\n");
    /* Append children of alpha just dequeued (if any) to end of queue. */
    if (alpha->child != NULL) {
      for (k=0; k<alpha->max_children; k++) {
        beta = alpha->child[k];
        if (beta != NULL) {
          enqueue (&q, beta);
          /* Print edge information into GraphViz file for the edge
             connecting alpha to beta. */
          fprintf (out_file, "%d->%d [label=\"%-8.4g\"];\n",
                   alpha->label,beta->label,beta->h);
        }
      }
    }
  }
  fprintf (out_file, "labelloc=\"t\";\n");
  fprintf (out_file, "label=\"%s\";\n", title);
  fprintf (out_file, "}\n");
  fclose (out_file);
  opts->tree_filename_num += 1;
  return;
}
