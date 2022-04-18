#include "test_pampac.h"


int
main () {
  options_struct opts = set_options (4,3);
  debug_print (3, &opts, __func__, 
               "opts.input_filename = %s\n", __func__, opts.input_filename);
  PTnode* root = make_tree (&opts);
  tree43 (root, &opts);
  opts.verbose = 5;
  assign_processes (root, &opts, 20);

  visualize_tree (root, &opts, "After prediction");
  prune_diverged_nodes (root, &opts);
  visualize_tree (root, &opts, "After correction");

  debug_print (0, &opts, __func__, "Cleaning up...\n");
  delete_tree (root, &opts);
  free_options (&opts);
  return 0;
}
