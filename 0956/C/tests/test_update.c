#include "test_pampac.h"

bool write_coordinates (int N, double *z) {
return true;
}

int
main () {
  /* Constructs long string of nodes to test advance_root_node. */
  options_struct opts = set_options (/*depth=*/5, /*width=*/1);
  opts.verbose = 5;
  PTnode* root = make_tree (&opts);
  /* Assign values to node fields. */
  root->child[0]->state = CONVERGED;
  root->child[0]->child[0]->state = CONVERGED;
  root->child[0]->child[0]->child[0]->state = CONVERGED;
  root->child[0]->child[0]->child[0]->child[0]->state = CONVERGING;
  assign_processes(root, &opts, 4);
  print_tree (root);

  visualize_tree(root, &opts, "Before advancing root");
  print_tree (root);
  advance_root_node (&root, &opts);

  visualize_tree(root, &opts, "After advancing root");
  print_tree (root);
  delete_tree(root, &opts);
  free_options (&opts);
  return 0;
}
