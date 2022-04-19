#include "pampac.h"
/**********************************************************************/
/* Loads data for initial point on curve from a file.                 */
/**********************************************************************/
bool
load_initial_coordinates (PTnode* node, options_struct* opts) {
  int N = opts->N_dim, k = 0;
  char *filename = opts->input_filename;
  double *z = node->z;
  FILE *input_file;

  debug_print (3, opts, __func__, "Reading 1st point to root node.\n");

  input_file = fopen (filename, "r");
  if (input_file == NULL) {
    debug_print (0, opts, __func__, "Error Opening File.\n");
    return (false);
  }
  /* Actually parse the input file. */
  for (k = 0; k < N; k++) {
    fscanf (input_file, "%lf", z);
    z++;
  }
  fclose (input_file);
  return (true);
}
