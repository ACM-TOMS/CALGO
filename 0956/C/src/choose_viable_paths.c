#include "pampac.h"
/**********************************************************************/
/* This function recursively traverses the tree identifying which     */
/* viable or valid paths to keep (deleting others in the process).    */
/**********************************************************************/
void
choose_viable_paths (PTnode *alpha, options_struct *opts) {
  int k, best_path_index;
  double valid_path_progress, viable_path_progress;
  PTnode *beta;

  /* At leaf nodes return immediately; otherwise, choose viable path */
  bool is_leaf = (count_children(alpha)==0);
  if (is_leaf)
    return;

  for (k=0; k<alpha->max_children; k++) {
    /* recursive call over child nodes */
    beta = alpha->child[k];
    if (beta!=NULL)
      choose_viable_paths (beta, opts);
  }

  bool has_valid_path = (alpha->valid_path_length>=0.0);
  bool has_viable_path = (alpha->viable_path_length>=0.0);
  best_path_index = -1;    /* Default: no subtress are deleted. */
  if (has_viable_path) {
    /* Default: if no valid path, use longest viable path. */
    best_path_index = alpha->viable_index;
    if (has_valid_path) {
      /* If viable path == valid path, no decision need be
         made; leave best_path_index as viable_path_index. */
      if (alpha->viable_index!=alpha->valid_index) {
        /* Principle decision: compare scaled rate of progress
        along viable and valid paths (if any exist). */
        viable_path_progress = (alpha->viable_path_length) /
                               (alpha->nu_viable+1.0);
        valid_path_progress  = (alpha->valid_path_length) /
                               (alpha->nu_valid+0.0);
        if (valid_path_progress >= viable_path_progress) {
          /* Choose valid path, delete viable path.
             Update fields of alpha accordingly. */
          best_path_index = alpha->valid_index;
          /* Note: valid is viable also by definition. */
          alpha->viable_index = alpha->valid_index;
          alpha->viable_path_length = alpha->valid_path_length;
        } else {
          /* Choose viable path, delete valid path.
             Update fields of alpha accordingly. */
          alpha->valid_index = -1;
          alpha->valid_path_length = NAN;
        }
      }
    }
  }
  if (best_path_index>=0) { /* A best path has been found */
    for (k=0; k<alpha->max_children; k++) {
      beta = alpha->child[k];
      /* Remove subtree rooted at node beta provided it has a
         viable path and it is *not* the best path. */
      if (beta!=NULL) {
        bool is_viable = (beta->viable_path_length>=0.0);
        if (is_viable && k!=best_path_index) {
          delete_tree (alpha->child[k], opts);
          alpha->child[k] = NULL;
        }
      }
    }
  }
  return;
}
