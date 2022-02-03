#include "pampac.h"
/**********************************************************************/
/* Blanks fields of options with default values (to make it easier in */
/* in principle to detect errors in which fields have not been        */
/* assigned). As such, pointers are initialised with NULL, integer-   */
/* valued fields are initialised with -1 and real-valued fields are   */
/* initialised with NAN.                                              */
/**********************************************************************/

bool validate_options (options_struct *opts) {
  bool is_valid = true;

  if (opts->N_dim<=0) {
    printf ("validate_options: N_DIM must be positive ");
    printf ("(N_DIM = %d)\n", opts->N_dim);
    is_valid = false;
  }

  if (isnan (opts->lambda_min)) {
    printf ("validate_options: LAMBDA_MIN not set ");
    printf ("(LAMBDA_MIN = %g)\n", opts->lambda_min);
    is_valid = false;
  }

  if (isnan (opts->lambda_max)) {
    printf ("validate_options: LAMBDA_MAX not set ");
    printf ("(LAMBDA_MAX = %g)\n", opts->lambda_max);
    is_valid = false;
  }

  if ((opts->lambda_index<0) || (opts->lambda_index>=opts->N_dim)) {
    printf ("validate_options: LAMBDA_INDEX ");
    printf ("must be between 0 and N_DIM-1 ");
    printf ("(LAMBDA_INDEX = %d)\n", opts->lambda_index);
    is_valid = false;
  }

  if ((opts->lambda_dir!=1) && (opts->lambda_dir!=-1)) {
    printf ("validate_options: LAMBDA_DIR ");
    printf ("must be +1 or -1 ");
    printf ("(LAMBDA_DIR = %d)\n", opts->lambda_dir);
    is_valid = false;
  }

  if (isnan (opts->delta_lambda)) {
    printf ("validate_options: DELTA_LAMBDA not set ");
    printf ("(DELTA_LAMBDA = %g)\n", opts->delta_lambda);
    is_valid = false;
  }

  if (isnan (opts->h_min) || opts->h_min<=0) {
    printf ("validate_options: H_MIN must be positive ");
    printf ("(H_MIN = %g)\n", opts->h_min);
    is_valid = false;
  }

  if (isnan (opts->h_max) || (opts->h_max<=0)) {
    printf ("validate_options: H_MAX must be positive ");
    printf ("(H_MAX = %g)\n", opts->h_max);
    is_valid = false;
  }

  if (isnan (opts->h_init)) {
    printf ("validate_options: H_INIT not set ");
    printf ("(H_INIT = %g)\n", opts->h_init);
    is_valid = false;
  }

  if (opts->max_iter<=0) {
    printf ("validate_options: MAX_ITER must be positive ");
    printf ("(MAX_ITER = %d)\n", opts->max_iter);
    is_valid = false;
  }

  if (isnan (opts->tol_residual) || (opts->tol_residual<=0)) {
    printf ("validate_options: TOL_RESIDUAL must be positive ");
    printf ("(TOL_RESIDUAL = %g)\n", opts->tol_residual);
    is_valid = false;
  }

  if (isnan (opts->mu) || (opts->mu <=0)) {
    printf ("validate_options: MU must be positive ");
    printf ("(MU = %g)\n", opts->mu);
    is_valid = false;
  }

  if (isnan (opts->gamma) || (opts->gamma<=0)) {
    printf ("validate_options: GAMMA must be positive ");
    printf ("(GAMMA = %g)\n", opts->gamma);
    is_valid = false;
  }

  if (opts->max_depth<=0) {
    printf ("validate_options: MAX_DEPTH must positive ");
    printf ("(MAX_DEPTH = %d)\n", opts->max_depth);
    is_valid = false;
  }

  if (opts->max_children<=0) {
    printf ("validate_options: MAX_CHILDREN must positive ");
    printf ("(MAX_CHILDREN = %d)\n", opts->max_children);
    is_valid = false;
  }

  if (opts->scale_factors==NULL) {
    printf("validate_options: SCALE_FACTORS not set\n");
    is_valid = false;
  } else {
    for (int k=0; k<opts->max_children; k++) {
      if (opts->scale_factors[k]<=0) {
        printf ("validate_options: SCALE_FACTORS ");
        printf ("must be positive ");
        printf ("(SCALE_FACTOR[%d] = %g)\n",
                k, opts->scale_factors[k]);
        is_valid = false;
      }
    }
  }

  if (opts->max_global_iter<=0) {
    printf("validate_options: MAX_GLOBAL_ITER must positive\n");
    printf ("(MAX_GLOBAL_ITER = %d)\n", opts->max_global_iter);
    is_valid = false;
  }

  if (opts->input_filename==NULL) {
    printf("validate_options: INPUT_FILENAME not set\n");
    is_valid = false;
  }

  if (opts->verbose<0) {
    printf("validate_options: VERBOSE not set\n");
    is_valid = false;
  }

  if ((opts->verbose>2) && (opts->tree_base_filename==NULL)) {
    printf("validate_options: TREE_BASE_FILENAME not set\n");
    is_valid = false;
  }

  return is_valid;
}
