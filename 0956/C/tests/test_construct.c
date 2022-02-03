#include "pampac.h"
extern options_struct set_options (int depth, int children);
extern PTnode* make_tree (options_struct* opts);

int
main () {
  options_struct opts = set_options (1,2);
  opts.verbose=5;
  PTnode* root = make_tree (&opts);
  opts.max_depth = 3;

  construct_predictor_nodes (root, &opts);
  assign_processes(root, &opts, 20);
  construct_predictor_nodes (root, &opts);
  assign_processes(root, &opts, 20);
  visualize_tree (root,&opts,"Before");

  delete_tree (root->child[0]->child[1], &opts);
  root->child[0]->child[1] = NULL;
  delete_tree (root->child[1]->child[0], &opts);
  root->child[1]->child[0] = NULL;
  delete_tree (root->child[1]->child[1], &opts);
  root->child[1]->child[1] = NULL;
  visualize_tree (root,&opts,"Middle");

  construct_predictor_nodes (root, &opts);
  assign_processes (root, &opts, 20);
  visualize_tree (root, &opts, "After");
  printf("Calling delete_tree...\n");
  delete_tree (root, &opts);
  free_options (&opts);
  return 0;
}
