#include <string.h>
#include "pampac.h"
/**********************************************************************/
/* Utility routine to translate strings parsed by parse_options into  */
/* values (and check validity) in the options data structure. A look- */
/* up table and an associated enumeration is introduced to make it    */
/* easier to modify (in the event of adding or removing options).     */
/**********************************************************************/
bool
assign_options (char *name, char* val, options_struct *opt) {
  enum  { OPT_N_DIM,
          OPT_LAMBDA_MIN,
          OPT_LAMBDA_MAX,
          OPT_LAMBDA_INDEX,
          OPT_DELTA_LAMBDA,
          OPT_H_MIN,
          OPT_H_MAX,
          OPT_H_INIT,
          OPT_MAX_ITER,
          OPT_TOL_RESIDUAL,
          OPT_GAMMA,
          OPT_MU,
          OPT_MAX_DEPTH,
          OPT_VERBOSE,
          OPT_MAX_GLOBAL_ITER,
          OPT_INPUT_FILENAME,
          OPT_TREE_BASE_FILENAME,
          OPT_SCALE_FACTOR
        };
  const char* table[] = {
    "N_DIM",
    "LAMBDA_MIN",
    "LAMBDA_MAX",
    "LAMBDA_INDEX",
    "DELTA_LAMBDA",
    "H_MIN",
    "H_MAX",
    "H_INIT",
    "MAX_ITER",
    "TOL_RESIDUAL",
    "GAMMA",
    "MU",
    "MAX_DEPTH",
    "VERBOSE",
    "MAX_GLOBAL_ITER",
    "INPUT_FILENAME",
    "TREE_BASE_FILENAME",
    "SCALE_FACTOR",
  };
  bool is_valid_name = true;
  const int N_options = sizeof table / sizeof *table;

  /* Linear search for name in the look-up table */
  int option_index = -1;
  for (int k=0; k != N_options; k++) {
    if ( strcmp( name, table[k]) == 0 ) {
      option_index = k;
      break;
    }
  }
  int n_chars;
  switch (option_index) {
  case OPT_N_DIM:
    opt->N_dim = atoi (val);
    break;
  case OPT_LAMBDA_MIN:
    opt->lambda_min = atof (val);
    break;
  case OPT_LAMBDA_MAX:
    opt->lambda_max = atof (val);
    break;
  case OPT_LAMBDA_INDEX:
    opt->lambda_index = atoi (val);
    break;
  case OPT_DELTA_LAMBDA:
    opt->delta_lambda = atof (val);
    break;
  case OPT_H_MIN:
    opt->h_min = atof (val);
    break;
  case OPT_H_MAX:
    opt->h_max = atof (val);
    break;
  case OPT_H_INIT:
    opt->h_init = atof (val);
    break;
  case OPT_MAX_ITER:
    opt->max_iter = atoi (val);
    break;
  case OPT_TOL_RESIDUAL:
    opt->tol_residual = atof (val);
    break;
  case OPT_GAMMA:
    opt->gamma = atof (val);
    break;
  case OPT_MU:
    opt->mu = atof (val);
    break;
  case OPT_MAX_DEPTH:
    opt->max_depth = atoi (val) + 1; /* User does not include root node in depth */
    break;
  case OPT_VERBOSE:
    opt->verbose = atoi (val);
    break;
  case OPT_MAX_GLOBAL_ITER:
    opt->max_global_iter = atoi (val);
    break;
  case OPT_INPUT_FILENAME:
    n_chars = strlen (val) + 1;
    opt->input_filename = malloc (n_chars * sizeof (char));
    memcpy (opt->input_filename, val, n_chars);
    break;
  case OPT_TREE_BASE_FILENAME:
    n_chars = strlen (val) + 1;
    opt->tree_base_filename = malloc (n_chars * sizeof (char));
    memcpy (opt->tree_base_filename, val, n_chars);
    break;
  case OPT_SCALE_FACTOR:
    if (opt->scale_factors==NULL)
      opt->scale_factors = malloc ( sizeof(double) );
    else
      opt->scale_factors = realloc ( opt->scale_factors,
                                     (1+opt->max_children) *
                                     sizeof (double) );
    double c = atof (val);
    opt->scale_factors[opt->max_children] = c;
    opt->max_children++;
    break;
  default:
    printf ("assign_options: Unable to match ");
    printf ("parameter name %s parsed.\n", name);
    is_valid_name = false;
  }

  return is_valid_name;
}
