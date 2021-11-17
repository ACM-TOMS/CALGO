#include "pampac.h"
/**********************************************************************/
/* Blanks fields of options with default values (to make it easier in */
/* in principle to detect uninitialised fields). As such, pointers    */
/* are initialised with NULL, integer-valued fields are initialised   */
/* with -1 and real-valued fields are initialised with NAN.           */
/**********************************************************************/
void initialize_options (options_struct *opts) {
  opts->N_dim = -1;
  opts->lambda_min = NAN;
  opts->lambda_max = NAN;
  opts->lambda_index = -1;
  opts->lambda_dir = 0;
  opts->delta_lambda = NAN;

  opts->h_min = NAN;
  opts->h_max = NAN;
  opts->h_init = NAN;

  opts->max_iter = -1;
  opts->tol_residual = NAN;
  opts->mu = NAN;
  opts->gamma = NAN;

  opts->max_depth = -1;
  opts->max_children = 0;
  opts->scale_factors = NULL;
  opts->max_global_iter = -1;

  opts->verbose = 0;
  opts->input_filename = NULL;
  opts->tree_base_filename = NULL;
  opts->tree_filename_num = 0;
  return;
}
