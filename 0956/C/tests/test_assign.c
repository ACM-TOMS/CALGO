#include "test_pampac.h"

int
main (int argc, char* argv[]) {
  options_struct opts = set_options (1,2);
  opts.verbose = 5;
  PTnode* root = make_tree (&opts);
  int nproc = atoi (argv[1]);
  opts.max_depth = 4;

  construct_predictor_nodes (root, &opts);
  assign_processes(root, &opts, nproc);
  assign_predictor_steps (root, &opts);
  visualize_tree(root, &opts, "After prediction");
  print_tree (root);
  printf("*********************************************************\n");

  construct_predictor_nodes (root, &opts);
  assign_processes (root, &opts, nproc);
  assign_predictor_steps (root, &opts);
  delete_tree (root->child[0]->child[1], &opts);
  root->child[0]->child[1] = NULL;
  delete_tree (root->child[1]->child[0], &opts);
  root->child[1]->child[0] = NULL;
  delete_tree (root->child[1]->child[1], &opts);
  root->child[1]->child[1] = NULL;
  printf("*********************************************************\n");

  construct_predictor_nodes (root, &opts);
  assign_processes(root, &opts, nproc);
  assign_predictor_steps (root, &opts);
  visualize_tree(root,&opts, "After prediction");
  print_tree (root);
  printf("Calling delete_tree...\n");
  delete_tree(root, &opts);
  return 0;
}
